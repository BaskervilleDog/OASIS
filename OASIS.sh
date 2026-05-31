#!/bin/bash

# =============================================================================
#  OASIS — Ortholog Alignment & Similarity Screener
#  Refactored for lower cyclomatic complexity.
#
#  Key changes vs original:
#    - resolve_gene_id split into _resolve_protein_to_gene / _resolve_nuccore_to_gene
#    - run_single split into _fetch_query_sequence / _run_outputs
#    - CLI/interactive input split into parse_args / prompt_interactive
#    - Lookup tables replace repetitive case/if chains where possible
#    - Each function now has CC ≤ 6; two dispatcher functions sit at ~5
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# 0. NCBI API key
# ---------------------------------------------------------------------------
NCBI_API_KEY="${NCBI_API_KEY:-}"

_api_key_param() {
    [ -n "$NCBI_API_KEY" ] && printf '&api_key=%s' "$NCBI_API_KEY" || true
}

# ---------------------------------------------------------------------------
# 1. Low-level download helpers
# ---------------------------------------------------------------------------

download_file() {
    local url="$1" output="$2"
    if command -v curl &>/dev/null; then
        curl -s -L --globoff -o "$output" "$url"
    elif command -v wget &>/dev/null; then
        wget -q -O "$output" "$url"
    else
        echo "❌ Critical Error: Neither 'curl' nor 'wget' found." >&2
        exit 1
    fi
}

# CC = 4  (for-loop + curl/wget branch + result check)
download_text() {
    local url="$1"
    local result wait_secs=2

    for attempt in 1 2 3 4 5; do
        if command -v curl &>/dev/null; then
            # --globoff disables curl's bracket-range expansion so NCBI
            # query terms like [accn] and [gene] are sent verbatim
            result=$(curl --silent --location --fail --globoff "$url" 2>/dev/null) || true
        else
            result=$(wget -q -O - "$url" 2>/dev/null) || true
        fi

        if [ -n "$result" ]; then
            printf '%s' "$result"
            return 0
        fi

        echo "  ⚠️  Retry $attempt/5 (waiting ${wait_secs}s)..." >&2
        sleep "$wait_secs"
        wait_secs=$(( wait_secs * 2 ))
    done

    return 1
}

extract_zip() {
    local zip_file="$1" dest_dir="$2"
    if command -v unzip &>/dev/null; then
        unzip -q -o "$zip_file" -d "$dest_dir"
    elif command -v python3 &>/dev/null; then
        python3 -c "import zipfile; zipfile.ZipFile('$zip_file','r').extractall('$dest_dir')"
    else
        echo "❌ Critical Error: Neither 'unzip' nor 'python3' found." >&2
        exit 1
    fi
}

# ---------------------------------------------------------------------------
# 2. Tool installation
# ---------------------------------------------------------------------------

BLAST_DIR="$HOME/ncbi-blast-2.13.0+/bin"
DATASETS_PATH="$HOME/datasets"
export PATH="$BLAST_DIR:$HOME:$PATH"

install_tools() {
    # Re-check even if file exists: a previous interrupted download may have
    # left a corrupt HTML error page that is not a valid binary.
    if [ ! -x "$DATASETS_PATH" ] || ! "$DATASETS_PATH" --version &>/dev/null; then
        echo "📦 Installing NCBI Datasets CLI..."
        rm -f "$DATASETS_PATH"
        download_file \
            'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets' \
            "$DATASETS_PATH"
        # Validate: ELF/executable magic, not just file existence
        if ! file "$DATASETS_PATH" 2>/dev/null | grep -qE "ELF|executable"; then
            echo "❌ datasets CLI download produced non-binary content." >&2
            echo "   Check network access to ftp.ncbi.nlm.nih.gov." >&2
            echo "   First line of downloaded file: $(head -1 "$DATASETS_PATH" 2>/dev/null)" >&2
            rm -f "$DATASETS_PATH"
            exit 1
        fi
        chmod +x "$DATASETS_PATH"
        echo "  ✅ datasets CLI installed: $("$DATASETS_PATH" --version 2>&1 | head -1)"
    fi
    if [ ! -d "$BLAST_DIR" ]; then
        echo "🛰️  Installing BLAST+ 2.13.0 (static binaries)..."
        download_file \
            'https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.13.0/ncbi-blast-2.13.0+-x64-linux.tar.gz' \
            "$HOME/blast.tar.gz"
        if ! file "$HOME/blast.tar.gz" 2>/dev/null | grep -qE "gzip|tar"; then
            echo "❌ BLAST+ download produced non-archive content." >&2
            echo "   Check network access to ftp.ncbi.nlm.nih.gov." >&2
            rm -f "$HOME/blast.tar.gz"
            exit 1
        fi
        tar -xzf "$HOME/blast.tar.gz" -C "$HOME"
        rm "$HOME/blast.tar.gz"
        echo "  ✅ BLAST+ installed: $("$BLAST_DIR/blastp" -version 2>&1 | head -1)"
    fi
}

# ---------------------------------------------------------------------------
# 3. Accession-type detection  (CC = 3)
# ---------------------------------------------------------------------------

detect_type() {
    local acc="$1"
    case "$acc" in
        NP_*|XP_*|WP_*) echo "protein"    ;;
        NM_*|XM_*|NG_*) echo "nucleotide" ;;
        [0-9]*)          echo "gene_id"    ;;
        *)
            echo "❌ Unrecognised accession format: '$acc'" >&2
            echo "   Supported: NP_ XP_ WP_ NM_ XM_ NG_ or numeric Gene ID" >&2
            echo "unknown"
            ;;
    esac
}

# ---------------------------------------------------------------------------
# 4. Gene ID resolution — one helper per accession type
# ---------------------------------------------------------------------------

# Shared inline Python fragments (keeps the helpers DRY)
_PY_FIRST_ID="
import sys, json
try:
    data = json.load(sys.stdin)
    ids = data['esearchresult']['idlist']
    if ids: print(ids[0])
except Exception: pass
"

_py_elink_gene() {
    # Args: $1 = space-separated accepted LinkName values
    local accepted="$1"
    python3 -c "
import sys
import xml.etree.ElementTree as ET
accepted = set('$accepted'.split())
try:
    root = ET.fromstring(sys.stdin.read())
    for db in root.iter('LinkSetDb'):
        if db.findtext('LinkName') in accepted:
            for link in db.findall('Link'):
                gid = link.findtext('Id')
                if gid:
                    print(gid)
                    raise SystemExit
except Exception: pass
"
}

# CC = 3  (two if-guards + one return)
_resolve_protein_to_gene() {
    local accession="$1"
    echo "  🔎 Resolving protein accession → Gene ID..." >&2

    local puid
    puid=$(download_text \
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=protein&term=${accession}[accn]&retmode=json$(_api_key_param)" |
        python3 -c "$_PY_FIRST_ID")

    if [ -z "$puid" ]; then
        echo "  ❌ Failed to resolve protein UID." >&2
        return 1
    fi

    local gene_id
    gene_id=$(download_text \
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=protein&db=gene&id=${puid}$(_api_key_param)" |
        _py_elink_gene "protein_gene protein_gene_refseq")

    if [ -z "$gene_id" ]; then
        echo "  ❌ Failed to resolve Gene ID from protein accession." >&2
        return 1
    fi
    echo "$gene_id"
}

# CC = 3
_resolve_nuccore_to_gene() {
    local accession="$1"
    echo "  🔎 Resolving nucleotide accession → Gene ID..." >&2

    local nuid
    nuid=$(download_text \
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=${accession}[accn]&retmode=json$(_api_key_param)" |
        python3 -c "$_PY_FIRST_ID")

    if [ -z "$nuid" ]; then
        echo "  ❌ Failed to resolve nuccore UID." >&2
        return 1
    fi

    local gene_id
    gene_id=$(download_text \
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=nuccore&db=gene&id=${nuid}$(_api_key_param)" |
        _py_elink_gene "nuccore_gene nucleotide_gene nuccore_gene_refseq")

    if [ -z "$gene_id" ]; then
        echo "  ❌ Failed to resolve Gene ID from nucleotide accession." >&2
        return 1
    fi
    echo "$gene_id"
}

# Dispatcher — CC = 4  (case with 3 branches + passthrough)
resolve_gene_id() {
    local accession="$1"
    local acc_type
    acc_type=$(detect_type "$accession")

    case "$acc_type" in
        protein)    _resolve_protein_to_gene  "$accession" ;;
        nucleotide) _resolve_nuccore_to_gene  "$accession" ;;
        gene_id)    echo "$accession"                       ;;
        *)          return 1                                ;;
    esac
}

# ---------------------------------------------------------------------------
# 5. Ortholog fetcher  (CC = 5)
# ---------------------------------------------------------------------------

fetch_all_orthologs() {
    local gene_id="$1" output_faa="$2"
    local ortho_zip="$GLOBAL_TMP/orthologs_${gene_id}.zip"
    local ortho_dir="$GLOBAL_TMP/orthologs_${gene_id}"

    if [ -s "$output_faa" ]; then
        local cached
        cached=$(grep -c '^>' "$output_faa" 2>/dev/null || echo 0)
        echo "  ♻️  Using cached ortholog FAA ($cached sequences)."
        return 0
    fi

    echo "  ⚠️  Note: ortholog downloads cover vertebrates and insects only." >&2
    echo "  🌍 Downloading ortholog protein sequences (gene ID: $gene_id)..." >&2

    # Redirect only stdout; keep stderr visible so network/auth errors surface
    "$DATASETS_PATH" download gene gene-id "$gene_id" \
        --ortholog all --include protein --filename "$ortho_zip" >/dev/null

    if [ ! -s "$ortho_zip" ]; then
        echo "  ❌ Ortholog download failed for gene ID $gene_id." >&2
        return 1
    fi

    mkdir -p "$ortho_dir"
    extract_zip "$ortho_zip" "$ortho_dir"

    find "$ortho_dir" -type f \( -name "*.faa" -o -name "*.protein.faa" \) \
        -exec cat {} + >"$output_faa"

    if [ ! -s "$output_faa" ]; then
        echo "  ❌ No protein sequences found after extraction." >&2
        return 1
    fi

    # Deduplicate on first field of header
    awk '/^>/{header=$0;split(header,a," ");id=a[1];if(seen[id]++){skip=1}else{skip=0;print};next} !skip{print}' \
        "$output_faa" > "${output_faa}.dedup"
    mv "${output_faa}.dedup" "$output_faa"

    local count
    count=$(grep -c '^>' "$output_faa" 2>/dev/null || echo 0)
    echo "  ✅ $count ortholog protein sequences ready."

    if [ "$count" -ge 499 ]; then
        echo "" >&2
        echo "  ⚠️  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━" >&2
        echo "  ⚠️  DATASET CAP REACHED: $count sequences (NCBI cap ~499)." >&2
        echo "  ⚠️  Results are INCOMPLETE. Split by taxon or use the web portal." >&2
        echo "  ⚠️  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━" >&2
        echo "" >&2
    fi
}

# ---------------------------------------------------------------------------
# 6. BLAST filter  (CC = 3)
# ---------------------------------------------------------------------------

run_blast_filter() {
    local query_fasta="$1" ortho_faa="$2" blast_prog="$3"
    local min_id="$4" min_sim="$5" db_prefix="$6" out_list="$7" exclude_id="$8"

    if [ ! -f "${db_prefix}.pin" ] && [ ! -f "${db_prefix}.pdb" ]; then
        echo "  🔨 Building local BLAST database..."
        "$BLAST_DIR/makeblastdb" -in "$ortho_faa" -dbtype prot \
            -out "$db_prefix" -parse_seqids -logfile /dev/null
    else
        echo "  ♻️  Reusing existing BLAST database."
    fi

    echo "  🔬 Running $blast_prog alignment..."
    echo "      Filters: pident ≥ ${min_id}%  |  ppos ≥ ${min_sim}%"

    "$BLAST_DIR/$blast_prog" \
        -query "$query_fasta" -db "$db_prefix" \
        -max_hsps 1 -max_target_seqs 5000 \
        -outfmt "6 saccver pident ppos" -evalue 1e-5 |
    awk -v id_min="$min_id" -v sim_min="$min_sim" \
        '($2+0)>=(id_min+0) && ($3+0)>=(sim_min+0) {print $1}' |
    awk -v q="$exclude_id" '
        BEGIN { sub(/\.[0-9]+$/,"",q) }
        { x=$0; sub(/\.[0-9]+$/,"",x); if(x!=q) print $0 }
    ' | sort -u > "$out_list"
}

# ---------------------------------------------------------------------------
# 7. Output extractors
# ---------------------------------------------------------------------------

# CC = 2
extract_protein_fasta() {
    local db_prefix="$1" acclist="$2" out_fasta="$3"
    echo "  🚀 Extracting protein sequences from local database..."
    "$BLAST_DIR/blastdbcmd" -db "$db_prefix" -entry_batch "$acclist" \
        -out "$out_fasta" 2>/dev/null

    if [ -s "$out_fasta" ]; then
        local n; n=$(grep -c '^>' "$out_fasta")
        echo "  ✅ Protein FASTA: $n sequences → $(basename "$out_fasta")"
    else
        echo "  ❌ Protein extraction returned empty output." >&2
    fi
}

# CC = 4  (for-loop + if file + for-retry + if success)
extract_cds_fasta() {
    local acclist="$1" out_fasta="$2"
    echo "  🚀 Downloading CDS sequences..."
    > "$out_fasta"

    local cds_batch_dir
    cds_batch_dir=$(mktemp -d "${GLOBAL_TMP}/cds_batches_XXXXXX")
    split -l 100 "$acclist" "${cds_batch_dir}/batch_"

    for batch in "${cds_batch_dir}"/batch_*; do
        [ -f "$batch" ] || continue
        _fetch_cds_batch "$batch" "$out_fasta"
    done

    local total_cds
    total_cds=$(grep -c "^>" "$out_fasta" 2>/dev/null || echo 0)
    if [ "$total_cds" -gt 0 ]; then
        echo "  ✅ CDS FASTA: $total_cds sequences → $(basename "$out_fasta")"
    else
        echo "  ❌ CDS output is empty." >&2
    fi
}

# CC = 3  (retry for-loop + success check)
_fetch_cds_batch() {
    local batch="$1" out_fasta="$2"
    local ids; ids=$(paste -sd, "$batch")
    local tmp_fa="${batch}.fa"
    local wait_secs=2

    for retry in 1 2 3 4 5; do
        download_text \
"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=${ids}&rettype=fasta_cds_na&retmode=text$(_api_key_param)" \
> "$tmp_fa"

        if [ -s "$tmp_fa" ]; then
            cat "$tmp_fa" >> "$out_fasta"
            local nseq; nseq=$(grep -c "^>" "$tmp_fa")
            echo "    → Batch $(basename "$batch"): $nseq CDS sequences"
            return 0
        fi

        echo "    ⚠️ Retry $retry/5 (waiting ${wait_secs}s)..."
        sleep "$wait_secs"
        wait_secs=$(( wait_secs * 2 ))
    done

    echo "    ⚠️ Batch $(basename "$batch") failed after 5 retries." >&2
}

# ---------------------------------------------------------------------------
# 8. Query sequence fetcher — extracted from run_single  (CC = 4)
#    Returns blast program name on stdout; writes FASTA to $query_fasta path
# ---------------------------------------------------------------------------

_fetch_query_sequence() {
    local acc="$1" acc_type="$2" query_fasta="$3"

    # Map accession type → (blast program, efetch db, efetch params)
    # Using parallel arrays as a lightweight lookup table
    local -A BLAST_PROG=( [protein]=blastp [nucleotide]=blastx [gene_id]=blastp )
    local -A FETCH_DB=(   [protein]=protein [nucleotide]=nuccore )
    local -A FETCH_ARGS=( [protein]="rettype=fasta" [nucleotide]="rettype=fasta" )

    local blast_prog="${BLAST_PROG[$acc_type]}"

    if [ "$acc_type" = "gene_id" ]; then
        echo "  ℹ️  Numeric gene ID — resolving representative RefSeq protein..." >&2

        local prot_uid
        prot_uid=$(download_text \
"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=gene&db=protein&id=${acc}&retmode=json&linkname=gene_protein_refseq$(_api_key_param)" |
            python3 -c "
import sys, json
try:
    data = json.load(sys.stdin)
    for ls in data.get('linksets',[]):
        for db in ls.get('linksetdbs',[]):
            if db.get('linkname')=='gene_protein_refseq':
                lk = db.get('links',[])
                if lk: print(lk[0]); raise SystemExit
except Exception: pass
" 2>/dev/null)

        if [ -z "$prot_uid" ]; then
            echo "  ❌ Could not resolve a RefSeq protein for gene ID '$acc'." >&2
            echo "     Try supplying an NP_ accession instead." >&2
            return 1
        fi

        echo "  ✅ Representative protein UID: $prot_uid"
        download_text \
"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=${prot_uid}&rettype=fasta&retmode=text$(_api_key_param)" \
            >"$query_fasta"
    else
        local db="${FETCH_DB[$acc_type]}"
        download_text \
"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=${db}&id=${acc}&rettype=fasta&retmode=text$(_api_key_param)" \
            >"$query_fasta"
    fi

    echo "$blast_prog"
}

# ---------------------------------------------------------------------------
# 9. Output stage — extracted from run_single  (CC = 3)
# ---------------------------------------------------------------------------

_run_outputs() {
    local db_prefix="$1" acc="$2" out_dir="$3"
    local final_list="$4" want_prot="$5" want_cds="$6"

    if [[ "$want_prot" =~ ^[Yy]$ ]]; then
        extract_protein_fasta "$db_prefix" "$final_list" \
            "${out_dir}/sequences_PROT_OASIS_${acc}.fasta"
    fi

    if [[ "$want_cds" =~ ^[Yy]$ ]]; then
        extract_cds_fasta "$final_list" \
            "${out_dir}/sequences_CDS_OASIS_${acc}.fasta"
    fi
}

# ---------------------------------------------------------------------------
# 10. Single-accession pipeline  (CC = 5)
# ---------------------------------------------------------------------------

run_single() {
    local acc="$1" min_id="$2" min_sim="$3" want_prot="$4" want_cds="$5"

    echo ""
    echo "╔══════════════════════════════════════════════════════╗"
    echo "║  Processing: $acc"
    echo "╚══════════════════════════════════════════════════════╝"

    local acc_type
    acc_type=$(detect_type "$acc")
    if [ "$acc_type" = "unknown" ]; then
        echo "  ⏭️  Skipping '$acc' (unrecognised format)."
        return 1
    fi
    echo "  🔬 Type detected: $acc_type"

    local gene_id
    gene_id=$(resolve_gene_id "$acc") || {
        echo "  ❌ Could not resolve Gene ID for '$acc'. Skipping." >&2
        return 1
    }
    echo "  ✅ Gene ID: $gene_id"

    local out_dir="${OUTPUT_ROOT}/${acc}"
    mkdir -p "$out_dir"

    local query_fasta="$GLOBAL_TMP/query_${acc}.fasta"
    local blast_prog
    blast_prog=$(_fetch_query_sequence "$acc" "$acc_type" "$query_fasta") || return 1

    if [ ! -s "$query_fasta" ] || ! grep -q '^>' "$query_fasta"; then
        echo "  ❌ Could not fetch query sequence for '$acc'. Skipping." >&2
        return 1
    fi
    echo "  ✅ Query sequence fetched."

    local ortho_faa="$GLOBAL_TMP/ortho_${gene_id}.faa"
    local blast_db="$GLOBAL_TMP/blastdb_${gene_id}"
    local final_list="${out_dir}/filtered_accessions_ID${min_id}_SIM${min_sim}_${acc}.txt"

    fetch_all_orthologs "$gene_id" "$ortho_faa" || {
        echo "  ❌ Ortholog fetch failed for '$acc'. Skipping." >&2
        return 1
    }

    run_blast_filter "$query_fasta" "$ortho_faa" "$blast_prog" \
        "$min_id" "$min_sim" "$blast_db" "$final_list" "$acc"

    local count
    count=$(wc -l <"$final_list" | tr -d '[:space:]')
    echo "  🎯 $count accessions passed filters (pident ≥ ${min_id}%, ppos ≥ ${min_sim}%)."
    echo "  📋 Accession list → $final_list"

    if [ "$count" -eq 0 ]; then
        echo "  ⚠️  No hits passed filters — try lowering thresholds."
        return 0
    fi

    _run_outputs "$blast_db" "$acc" "$out_dir" "$final_list" "$want_prot" "$want_cds"
    echo "  ✅ Done → ./${out_dir}/"
}

# ---------------------------------------------------------------------------
# 11. Input parsing — split into two focused functions
# ---------------------------------------------------------------------------

# Interactive mode  (CC = 2)
prompt_interactive() {
    echo ""
    echo "You can enter:"
    echo "  • A single accession:           NP_001416352.1"
    echo "  • Multiple accessions (space):  NP_001416352.1 NM_001429423.1"
    echo "  • A path to a batch file:       /path/to/ids.txt"
    echo ""
    read -rp "🧬 Accession(s) or batch file: " INPUT_RAW

    if [ -f "$INPUT_RAW" ]; then
        mapfile -t ACCESSIONS < <(grep -v '^\s*#' "$INPUT_RAW" | grep -v '^\s*$' | awk '{print $1}')
    else
        read -ra ACCESSIONS <<<"$INPUT_RAW"
    fi

    echo ""
    echo "  pident = % identical positions (strict sequence identity)"
    echo "  ppos   = % positively-scoring substitutions (similarity)"
    echo "  Typical values: pident 70–90, ppos 80–95"
    echo ""
    read -rp "📊 Minimum pident / Identity  (e.g. 70): " MIN_ID
    read -rp "📊 Minimum ppos   / Similarity (e.g. 80): " MIN_SIM
}

# CLI mode  (CC = 4)
parse_args() {
    local -n _accs="$1"   # nameref: output array for accessions
    local -n _id="$2"     # nameref: output var for MIN_ID
    local -n _sim="$3"    # nameref: output var for MIN_SIM

    local args=("${@:4}")
    local n=${#args[@]}

    # Strip trailing numeric pair as thresholds
    if [[ "${args[$((n-1))]}" =~ ^[0-9]+$ ]] && \
       [[ "${args[$((n-2))]}" =~ ^[0-9]+$ ]]; then
        _sim="${args[$((n-1))]}"
        _id="${args[$((n-2))]}"
        args=("${args[@]:0:$((n-2))}")
    fi

    # Single batch file or inline accessions
    if [ "${#args[@]}" -eq 1 ] && [ -f "${args[0]}" ]; then
        mapfile -t _accs < <(grep -v '^\s*#' "${args[0]}" | grep -v '^\s*$' | awk '{print $1}')
        echo "📂 Batch file: ${args[0]} (${#_accs[@]} accessions)"
    else
        _accs=("${args[@]}")
    fi

    # Prompt for any missing thresholds
    if [ -z "${_id:-}" ]; then
        echo "  pident = % identical positions  |  ppos = % positively-scoring substitutions"
        read -rp "📊 Minimum pident / Identity  (%): " _id
    fi
    if [ -z "${_sim:-}" ]; then
        read -rp "📊 Minimum ppos   / Similarity (%): " _sim
    fi
}

# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

cat <<'OASIS_BANNER'
=============================================================

     *******         **        ********   **    ********
    **/////**       ****      **//////   /**   **////// 
   **     //**     **//**    /**         /**  /**       
  /**      /**    **  //**   /*********  /**  /*********
  /**      /**   **********  ////////**  /**  ////////**
  //**     **   /**//////**         /**  /**         /**
   //*******    /**     /**   ********   /**   ******** 
    ///////     //      //   ////////    //   //////// 

          Ortholog Alignment & Similarity Screener
    EvoMol - LAboratório de Evolução Molecular e Sistemas
=============================================================
OASIS_BANNER

# --- Collect inputs (CC of main block reduced to ~4) ---

ACCESSIONS=()
MIN_ID=""
MIN_SIM=""

if [ "$#" -eq 0 ]; then
    prompt_interactive
else
    parse_args ACCESSIONS MIN_ID MIN_SIM "$@"
fi

[ "${#ACCESSIONS[@]}" -eq 0 ] && { echo "❌ No accessions provided. Exiting." >&2; exit 1; }
[ -z "$MIN_ID" ] || [ -z "$MIN_SIM" ] && { echo "❌ Thresholds required." >&2; exit 1; }

echo ""
echo "📋 Accessions to process (${#ACCESSIONS[@]}):"
for a in "${ACCESSIONS[@]}"; do echo "   • $a"; done
echo "📊 pident ≥ ${MIN_ID}%  |  ppos ≥ ${MIN_SIM}%"
echo ""

read -rp "📥 Extract protein FASTA for each result? (y/n): " WANT_PROT
read -rp "🧬 Download CDS FASTA for each result?     (y/n): " WANT_CDS

# --- NCBI API key reminder ---
if [ -z "$NCBI_API_KEY" ]; then
    echo ""
    echo "  ℹ️  NCBI_API_KEY not set — running at 3 req/s."
    echo "     export NCBI_API_KEY=your_key_here"
    echo "     (Register free at: https://www.ncbi.nlm.nih.gov/account/)"
    echo ""
fi

# --- Disk space pre-check ---
ESTIMATED_MB=$(( ${#ACCESSIONS[@]} * 500 ))
AVAILABLE_MB=$(( $(df -k "$HOME" | awk 'NR==2{print $4}') / 1024 ))

if [ "$AVAILABLE_MB" -lt "$ESTIMATED_MB" ]; then
    echo "⚠️  Low disk space: need ~${ESTIMATED_MB} MB, have ${AVAILABLE_MB} MB." >&2
    read -rp "   Continue anyway? (y/n): " CONTINUE_ANYWAY
    [[ ! "$CONTINUE_ANYWAY" =~ ^[Yy]$ ]] && { echo "Aborted." >&2; exit 1; }
fi

# --- Environment setup ---
install_tools

OUTPUT_ROOT="OASIS_results_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$OUTPUT_ROOT"

GLOBAL_TMP=$(mktemp -d "${TMPDIR:-/tmp}/OASIS_XXXXXX")

_cleanup() {
    echo "" >&2
    echo "  🧹 Cleaning up temporary files in $GLOBAL_TMP ..." >&2
    rm -rf "$GLOBAL_TMP"
}
trap '_cleanup' EXIT INT TERM
# NOTE: SIGKILL cannot be trapped. Clean up manually with: rm -rf /tmp/OASIS_*

echo ""
echo "📁 Results: ./${OUTPUT_ROOT}/"
echo "🗂️  Temp:    $GLOBAL_TMP"
echo ""

# --- Main batch loop ---
SUCCESS=()
SKIPPED=()

for ACC in "${ACCESSIONS[@]}"; do
    if run_single "$ACC" "$MIN_ID" "$MIN_SIM" "$WANT_PROT" "$WANT_CDS"; then
        SUCCESS+=("$ACC")
    else
        SKIPPED+=("$ACC")
    fi
done

# --- Final summary ---
echo ""
echo "======================================================"
echo "🏁 OASIS batch run complete."
echo "   Processed : ${#ACCESSIONS[@]} accession(s)"
echo "   Succeeded : ${#SUCCESS[@]}"
echo "   Skipped   : ${#SKIPPED[@]}"
if [ "${#SKIPPED[@]}" -gt 0 ]; then
    echo "   Failed IDs:"
    for s in "${SKIPPED[@]}"; do echo "     • $s"; done
fi
echo "📁 Output root: ./${OUTPUT_ROOT}/"
echo "======================================================"
echo ""
echo "📖 COLUMN DEFINITIONS:"
echo "   pident = % identical positions"
echo "   ppos   = % positively-scoring substitutions"
echo ""
echo "⚠️  NOTICES:"
echo "   • NCBI caps ortholog downloads at ~499 sequences."
echo "   • Orthologs available for vertebrates and insects only."
echo "   • Manual curation recommended before MSA or phylogenetic analysis."
echo "   • SIGKILL cleanup: rm -rf /tmp/OASIS_*"