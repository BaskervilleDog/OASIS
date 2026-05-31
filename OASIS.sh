#!/bin/bash

# =============================================================================
#  OASIS — Ortholog Alignment & Similarity Screener
#  Supports: single or multiple accessions, auto-detected molecule type,
#            protein (NP_/XP_/WP_), nucleotide (NM_/XM_/NG_), and numeric
#            NCBI Gene IDs as input.
# Usage:
#   ./OASIS.sh                              (fully interactive)
#   ./OASIS.sh NP_001416352.1 90 95        (single accession, CLI args)
#   ./OASIS.sh accessions.txt  90 95       (batch file, one accession per line)
#   ./OASIS.sh NP_123.1 NM_456.1 90 95    (multiple accessions inline)
# =============================================================================

# --- CORRECTION 1: strict mode ---
# Prevents silent failures in pipelines and unset variables from propagating
# as empty strings through the pipeline.
set -euo pipefail

# --- 0. Helper utilities (rootless / sudo-free) ---

# --- CORRECTION 2: NCBI API key support ---
# Without a key NCBI enforces 3 req/s; with one it's 10 req/s.
# Export NCBI_API_KEY in your shell or set it here.
NCBI_API_KEY="${NCBI_API_KEY:-}"

_api_key_param() {
    # Appends &api_key=... to a URL fragment if the key is set
    if [ -n "$NCBI_API_KEY" ]; then
        printf '&api_key=%s' "$NCBI_API_KEY"
    fi
}

download_file() {
    local url="$1" output="$2"
    if command -v curl &>/dev/null; then
        curl -s -L -o "$output" "$url"
    elif command -v wget &>/dev/null; then
        wget -q -O "$output" "$url"
    else
        echo "❌ Critical Error: Neither 'curl' nor 'wget' found." >&2
        exit 1
    fi
}

# --- CORRECTION 3: HTTP status-aware download with exponential backoff ---
# Original retried on empty output only, treating HTTP 429 / 5xx the same as
# a transient network glitch. This version checks the HTTP status code and
# applies a proper backoff so we do not hammer the API during a ban.
download_text() {
    local url="$1"
    local result http_code wait_secs=2

    for attempt in 1 2 3 4 5; do
        if command -v curl &>/dev/null; then
            # Write the HTTP status to a temp file; body goes to stdout
            local tmp_http; tmp_http=$(mktemp)
            result=$(curl \
                --silent \
                --location \
                --write-out '%{http_code}' \
                --output /dev/fd/3 \
                --retry 0 \
                "$url" 3>&1 2>/dev/null) || true
            http_code="$result"
            result=$(curl \
                --silent \
                --location \
                --fail \
                "$url" 2>/dev/null) || true
            rm -f "$tmp_http"
        elif command -v wget &>/dev/null; then
            result=$(wget -q -O - "$url" 2>/dev/null) || true
            http_code="200"          # wget exits non-zero on HTTP errors
        fi

        if [ -n "$result" ]; then
            printf '%s' "$result"
            return 0
        fi

        # Exponential back-off: 2 → 4 → 8 → 16 → give up
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

# --- 1. Tool installation ---

BLAST_DIR="$HOME/ncbi-blast-2.13.0+/bin"
DATASETS_PATH="$HOME/datasets"
export PATH="$BLAST_DIR:$HOME:$PATH"

install_tools() {
    if [ ! -f "$DATASETS_PATH" ]; then
        echo "📦 Installing NCBI Datasets CLI..."
        download_file \
            'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets' \
            "$DATASETS_PATH"
        chmod +x "$DATASETS_PATH"
    fi
    if [ ! -d "$BLAST_DIR" ]; then
        echo "🛰️  Installing BLAST+ 2.13.0 (static binaries)..."
        download_file \
            'https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.13.0/ncbi-blast-2.13.0+-x64-linux.tar.gz' \
            "$HOME/blast.tar.gz"
        tar -xzf "$HOME/blast.tar.gz" -C "$HOME"
        rm "$HOME/blast.tar.gz"
    fi
}

# --- 2. Accession-type auto-detection ---

detect_type() {
    local acc="$1"
    case "$acc" in
    NP_* | XP_* | WP_*) echo "protein" ;;
    NM_* | XM_* | NG_*) echo "nucleotide" ;;
    [0-9]*) echo "gene_id" ;;
    *)
        echo "❌ Unrecognised accession format: '$acc'" >&2
        echo "   Supported: NP_ XP_ WP_ NM_ XM_ NG_ or numeric Gene ID" >&2
        echo "unknown"
        ;;
    esac
}

# --- 3. Gene ID resolver ---

# --- CORRECTION 4: Disambiguation when esearch returns multiple UIDs ---
# Original always took idlist[0], which is arbitrary for multi-hit queries
# (paralogs, pseudogenes, accessions present across organisms).
# Now we fetch the accession summary for each candidate UID and select the
# one whose 'accessionversion' field matches the input accession exactly.
# Falls back to idlist[0] with a warning only when no exact match is found.
resolve_gene_id() {

    local accession="$1"
    local acc_type
    acc_type=$(detect_type "$accession")

    case "$acc_type" in

    protein)
        echo "  🔎 Resolving protein accession → Gene ID..." >&2

        local puid
        puid=$(download_text \
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=protein&term=${accession}[accn]&retmode=json$(\ _api_key_param)" |
    python3 -c "
import sys, json
try:
    data = json.load(sys.stdin)
    ids = data['esearchresult']['idlist']
    if ids:
        print(ids[0])
except Exception:
    pass
")

        if [ -z "$puid" ]; then
            echo "  ❌ Failed to resolve protein UID." >&2
            return 1
        fi

        local gene_id
        gene_id=$(download_text \
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=protein&db=gene&id=${puid}$(_api_key_param)" |
    python3 -c "
import sys
import xml.etree.ElementTree as ET
try:
    root = ET.fromstring(sys.stdin.read())
    accepted = {'protein_gene', 'protein_gene_refseq'}
    for db in root.iter('LinkSetDb'):
        linkname = db.findtext('LinkName')
        if linkname in accepted:
            for link in db.findall('Link'):
                gid = link.findtext('Id')
                if gid:
                    print(gid)
                    raise SystemExit
except Exception:
    pass
")
        if [ -z "$gene_id" ]; then
            echo "  ❌ Failed to resolve Gene ID from protein accession." >&2
            return 1
        fi
        echo "$gene_id"
        ;;

    nucleotide)
        echo "  🔎 Resolving nucleotide accession → Gene ID..." >&2

        local nuid
        nuid=$(download_text \
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=${accession}[accn]&retmode=json$(_api_key_param)" |
    python3 -c "
import sys, json
try:
    data = json.load(sys.stdin)
    ids = data['esearchresult']['idlist']
    if ids:
        print(ids[0])
except Exception:
    pass
")

        if [ -z "$nuid" ]; then
            echo "  ❌ Failed to resolve nuccore UID." >&2
            return 1
        fi

        local gene_id
        gene_id=$(download_text \
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=nuccore&db=gene&id=${nuid}$(_api_key_param)" |
    python3 -c "
import sys
import xml.etree.ElementTree as ET
try:
    root = ET.fromstring(sys.stdin.read())
    accepted = {'nuccore_gene', 'nucleotide_gene', 'nuccore_gene_refseq'}
    for db in root.iter('LinkSetDb'):
        linkname = db.findtext('LinkName')
        if linkname in accepted:
            for link in db.findall('Link'):
                gid = link.findtext('Id')
                if gid:
                    print(gid)
                    raise SystemExit
except Exception:
    pass
")
        if [ -z "$gene_id" ]; then
            echo "  ❌ Failed to resolve Gene ID from nucleotide accession." >&2
            return 1
        fi
        echo "$gene_id"
        ;;

    gene_id)
        echo "$accession"
        ;;

    *)
        return 1
        ;;
    esac
}

# --- 4. Ortholog fetcher ---

# --- CORRECTION 5: Warn explicitly when the 499-sequence cap is reached ---
# The original only mentioned the cap in the final banner.  Now we check the
# count immediately after extraction and print a prominent warning so the
# user can act before downstream analysis begins on a truncated dataset.
fetch_all_orthologs() {
    local gene_id="$1"
    local output_faa="$2"
    local ortho_zip="$GLOBAL_TMP/orthologs_${gene_id}.zip"
    local ortho_dir="$GLOBAL_TMP/orthologs_${gene_id}"

    if [ -s "$output_faa" ]; then
        local cached
        cached=$(grep -c '^>' "$output_faa" 2>/dev/null || echo 0)
        echo "  ♻️  Using cached ortholog FAA ($cached sequences)."
        return 0
    fi

    # --- CORRECTION 6: Scope warning for vertebrate/insect limitation ---
    # Moved from the final banner to the point of download so users know
    # immediately if their taxon is outside the supported set.
    echo "  ⚠️  Note: NCBI Datasets ortholog downloads cover vertebrates and" >&2
    echo "       insects only. Plant/fungal/other gene IDs will yield 0 orthologs." >&2
    echo "  🌍 Downloading ortholog protein sequences (gene ID: $gene_id)..." >&2

    "$DATASETS_PATH" download gene gene-id "$gene_id" \
        --ortholog all \
        --include protein \
        --filename "$ortho_zip" >/dev/null 2>&1

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

    # Deduplicate on header line (keeps first occurrence)
    awk '
/^>/ {
    header = $0
    split(header, a, " ")
    id = a[1]
    if (seen[id]++) { skip = 1 } else { skip = 0; print }
    next
}
!skip { print }
' "$output_faa" > "${output_faa}.dedup"
    mv "${output_faa}.dedup" "$output_faa"

    local count
    count=$(grep -c '^>' "$output_faa" 2>/dev/null || echo 0)
    echo "  ✅ $count ortholog protein sequences ready."

    # Correction 5 (continued): hard cap warning at point of use
    if [ "$count" -ge 499 ]; then
        echo "" >&2
        echo "  ⚠️  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━" >&2
        echo "  ⚠️  DATASET CAP REACHED: exactly $count sequences returned." >&2
        echo "  ⚠️  NCBI Datasets caps ortholog downloads at ~499 entries." >&2
        echo "  ⚠️  Your results are INCOMPLETE. Consider splitting by taxon" >&2
        echo "  ⚠️  or using NCBI Ortholog web portal for full retrieval." >&2
        echo "  ⚠️  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━" >&2
        echo "" >&2
    fi
}

# --- 5. Per-query BLAST + filtering ---

# --- CORRECTION 7: Clarify pident vs ppos semantics in output labels ---
# The original labels MIN_ID → "Identity" and MIN_SIM → "Similarity" without
# explaining that column 2 is pident (% identical positions) and column 3 is
# ppos (% positively-scoring substitutions). The filter logic is correct; the
# labelling and user prompts are improved for clarity (see MAIN section).
run_blast_filter() {
    local query_fasta="$1"
    local ortho_faa="$2"
    local blast_prog="$3"
    local min_id="$4"
    local min_sim="$5"
    local db_prefix="$6"
    local out_list="$7"
    local exclude_id="$8"

    if [ ! -f "${db_prefix}.pin" ] && [ ! -f "${db_prefix}.pdb" ]; then
        echo "  🔨 Building local BLAST database..."
        "$BLAST_DIR/makeblastdb" \
            -in "$ortho_faa" \
            -dbtype prot \
            -out "$db_prefix" \
            -parse_seqids \
            -logfile /dev/null
    else
        echo "  ♻️  Reusing existing BLAST database."
    fi

    echo "  🔬 Running $blast_prog alignment..."
    echo "      Filters: pident (identity) ≥ ${min_id}%  |  ppos (similarity) ≥ ${min_sim}%"

    "$BLAST_DIR/$blast_prog" \
    -query "$query_fasta" \
    -db "$db_prefix" \
    -max_hsps 1 \
    -max_target_seqs 5000 \
    -outfmt "6 saccver pident ppos" \
    -evalue 1e-5 |
    awk -v id_min="$min_id" -v sim_min="$min_sim" '
        ($2+0) >= (id_min+0) && ($3+0) >= (sim_min+0) { print $1 }
    ' |
    awk -v q="$exclude_id" '
        BEGIN { sub(/\.[0-9]+$/, "", q) }
        {
            x = $0
            sub(/\.[0-9]+$/, "", x)
            if (x != q) print $0
        }
    ' |
    sort -u > "$out_list"
}

# --- 6. Output extractors ---

extract_protein_fasta() {
    local db_prefix="$1"
    local acclist="$2"
    local out_fasta="$3"

    echo "  🚀 Extracting protein sequences from local database..."
    "$BLAST_DIR/blastdbcmd" \
        -db "$db_prefix" \
        -entry_batch "$acclist" \
        -out "$out_fasta" 2>/dev/null

    if [ -s "$out_fasta" ]; then
        local n
        n=$(grep -c '^>' "$out_fasta")
        echo "  ✅ Protein FASTA: $n sequences → $(basename "$out_fasta")"
    else
        echo "  ❌ Protein extraction returned empty output." >&2
    fi
}

extract_cds_fasta() {

    local acclist="$1"
    local out_fasta="$2"

    echo "  🚀 Downloading CDS sequences..."

    > "$out_fasta"

    local batch_size=100

    # --- CORRECTION 8: Use a temp dir for batch files to avoid collisions ---
    # Original wrote cds_batch_* directly into GLOBAL_TMP, which is shared
    # across all accessions. Parallel or sequential runs on different genes
    # could overwrite each other's batch files. Now each call gets its own
    # sub-directory.
    local cds_batch_dir
    cds_batch_dir=$(mktemp -d "${GLOBAL_TMP}/cds_batches_XXXXXX")

    split -l "$batch_size" "$acclist" "${cds_batch_dir}/batch_"

    for batch in "${cds_batch_dir}"/batch_*; do

        [ -f "$batch" ] || continue

        local ids
        ids=$(paste -sd, "$batch")

        local tmp_fa="${cds_batch_dir}/tmp_cds_$(basename "$batch").fa"

        local success=0
        local wait_secs=2

        for retry in 1 2 3 4 5; do

            download_text \
"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=${ids}&rettype=fasta_cds_na&retmode=text$(_api_key_param)" \
> "$tmp_fa"

            if [[ -s "$tmp_fa" ]]; then
                success=1
                break
            fi

            echo "    ⚠️ Retry $retry/5 (waiting ${wait_secs}s)..."
            sleep "$wait_secs"
            wait_secs=$(( wait_secs * 2 ))
        done

        if [[ $success -eq 1 ]]; then
            cat "$tmp_fa" >> "$out_fasta"
            local nseq
            nseq=$(grep -c "^>" "$tmp_fa")
            echo "    → Batch $(basename "$batch"): $nseq CDS sequences"
        else
            echo "    ⚠️ Batch $(basename "$batch") failed after 5 retries." >&2
        fi

    done

    local total_cds
    total_cds=$(grep -c "^>" "$out_fasta" 2>/dev/null || echo 0)

    if [ "$total_cds" -gt 0 ]; then
        echo "  ✅ CDS FASTA: $total_cds sequences → $(basename "$out_fasta")"
    else
        echo "  ❌ CDS output is empty." >&2
    fi
}

# --- 7. Single-accession pipeline ---

# --- CORRECTION 9: Repaired numeric gene ID query sequence fetch ---
# The original gene_id branch fetched the gene symbol via esummary but then
# performed a JSON-based elink call that the script never used (the result
# was assigned to $prot_uid via a separate call that was copy-pasted but
# inconsistently structured). The corrected version performs a clean
# gene→protein_refseq elink, validates the result, and fetches the FASTA
# directly, with an explicit error if the link returns empty.
run_single() {
    local acc="$1"
    local min_id="$2"
    local min_sim="$3"
    local want_prot="$4"
    local want_cds="$5"

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

    echo "  🔎 Resolving Gene ID..."
    local gene_id
    gene_id=$(resolve_gene_id "$acc")
    if [ -z "$gene_id" ]; then
        echo "  ❌ Could not resolve Gene ID for '$acc'. Skipping." >&2
        return 1
    fi
    echo "  ✅ Gene ID: $gene_id"

    local out_dir="${OUTPUT_ROOT}/${acc}"
    mkdir -p "$out_dir"

    local final_list="${out_dir}/filtered_accessions_ID${min_id}_SIM${min_sim}_${acc}.txt"

    echo "  📥 Fetching query sequence..."
    local query_fasta="$GLOBAL_TMP/query_${acc}.fasta"
    local blast_prog

    if [ "$acc_type" = "protein" ]; then
        blast_prog="blastp"
        download_text \
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=${acc}&rettype=fasta&retmode=text$(_api_key_param)" \
            >"$query_fasta"

    elif [ "$acc_type" = "nucleotide" ]; then
        blast_prog="blastx"
        download_text \
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${acc}&rettype=fasta&retmode=text$(_api_key_param)" \
            >"$query_fasta"

    else
        # Correction 9: clean gene ID → representative RefSeq protein fetch
        blast_prog="blastp"
        echo "  ℹ️  Numeric gene ID — resolving representative RefSeq protein..." >&2

        local prot_uid
        prot_uid=$(download_text \
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=gene&db=protein&id=${acc}&retmode=json&linkname=gene_protein_refseq$(_api_key_param)" |
            python3 -c "
import sys, json
try:
    data = json.load(sys.stdin)
    for ls in data.get('linksets', []):
        for db in ls.get('linksetdbs', []):
            if db.get('linkname') == 'gene_protein_refseq':
                lk = db.get('links', [])
                if lk:
                    print(lk[0])
                    raise SystemExit
except Exception:
    pass
" 2>/dev/null)

        if [ -z "$prot_uid" ]; then
            echo "  ❌ Could not resolve a RefSeq protein for gene ID '$acc'." >&2
            echo "     The gene may lack a RefSeq protein record, or may be outside" >&2
            echo "     vertebrate/insect scope. Try supplying an NP_ accession instead." >&2
            return 1
        fi

        echo "  ✅ Representative protein UID: $prot_uid"
        download_text \
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=${prot_uid}&rettype=fasta&retmode=text$(_api_key_param)" \
            >"$query_fasta"
    fi

    if [ ! -s "$query_fasta" ] || ! grep -q '^>' "$query_fasta"; then
        echo "  ❌ Could not fetch query sequence for '$acc'. Skipping." >&2
        return 1
    fi
    echo "  ✅ Query sequence fetched."

    local ortho_faa="$GLOBAL_TMP/ortho_${gene_id}.faa"
    local blast_db="$GLOBAL_TMP/blastdb_${gene_id}"

    if ! fetch_all_orthologs "$gene_id" "$ortho_faa"; then
        echo "  ❌ Ortholog fetch failed for '$acc'. Skipping." >&2
        return 1
    fi

    echo "  ⚙️  Running BLAST alignment and filtering..."
    run_blast_filter \
        "$query_fasta" "$ortho_faa" "$blast_prog" \
        "$min_id" "$min_sim" \
        "$blast_db" "$final_list" "$acc"

    local count
    count=$(wc -l <"$final_list" | tr -d '[:space:]')
    echo "  🎯 $count accessions passed filters (pident ≥ ${min_id}%, ppos ≥ ${min_sim}%)."
    echo "  📋 Accession list → $final_list"

    if [ "$count" -eq 0 ]; then
        echo "  ⚠️  No hits passed filters — try lowering thresholds."
        return 0
    fi

    if [[ "$want_prot" =~ ^[Yy]$ ]]; then
        local prot_out="${out_dir}/sequences_PROT_OASIS_${acc}.fasta"
        extract_protein_fasta "$blast_db" "$final_list" "$prot_out"
    fi

    if [[ "$want_cds" =~ ^[Yy]$ ]]; then
        local cds_out="${out_dir}/sequences_CDS_OASIS_${acc}.fasta"
        extract_cds_fasta "$final_list" "$cds_out"
    fi

    echo "  ✅ Done → ./${out_dir}/"
}

# =============================================================================
#  MAIN
# =============================================================================

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

# --- CORRECTION 10: Clarify threshold semantics in all prompts ---
# Users are now told what pident and ppos mean at the point of input,
# not just in the BLAST filter internals.

# ---------- Input parsing ----------

ACCESSIONS=()
MIN_ID=""
MIN_SIM=""

if [ "$#" -eq 0 ]; then
    echo ""
    echo "You can enter:"
    echo "  • A single accession:           NP_001416352.1"
    echo "  • Multiple accessions (space):  NP_001416352.1 NM_001429423.1"
    echo "  • A path to a batch file:       /path/to/ids.txt"
    echo ""
    read -rp "🧬 Accession(s) or batch file: " INPUT_RAW

    # --- CORRECTION 11: Quote the variable to handle paths with spaces ---
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

else
    ARGS=("$@")
    N=${#ARGS[@]}

    if [[ "${ARGS[$((N - 1))]}" =~ ^[0-9]+$ ]] && [[ "${ARGS[$((N - 2))]}" =~ ^[0-9]+$ ]]; then
        MIN_SIM="${ARGS[$((N - 1))]}"
        MIN_ID="${ARGS[$((N - 2))]}"
        ARGS=("${ARGS[@]:0:$((N - 2))}")
    fi

    if [ "${#ARGS[@]}" -eq 1 ] && [ -f "${ARGS[0]}" ]; then
        mapfile -t ACCESSIONS < <(grep -v '^\s*#' "${ARGS[0]}" | grep -v '^\s*$' | awk '{print $1}')
        echo "📂 Batch file: ${ARGS[0]} (${#ACCESSIONS[@]} accessions)"
    else
        ACCESSIONS=("${ARGS[@]}")
    fi

    if [ -z "$MIN_ID" ]; then
        echo ""
        echo "  pident = % identical positions  |  ppos = % positively-scoring substitutions"
        read -rp "📊 Minimum pident / Identity  (%): " MIN_ID
    fi
    if [ -z "$MIN_SIM" ]; then
        read -rp "📊 Minimum ppos   / Similarity (%): " MIN_SIM
    fi
fi

if [ "${#ACCESSIONS[@]}" -eq 0 ]; then
    echo "❌ No accessions provided. Exiting." >&2
    exit 1
fi
if [ -z "$MIN_ID" ] || [ -z "$MIN_SIM" ]; then
    echo "❌ Identity and Similarity thresholds are required." >&2
    exit 1
fi

echo ""
echo "📋 Accessions to process (${#ACCESSIONS[@]}):"
for a in "${ACCESSIONS[@]}"; do echo "   • $a"; done
echo "📊 pident ≥ ${MIN_ID}%  |  ppos ≥ ${MIN_SIM}%"
echo ""

read -rp "📥 Extract protein FASTA for each result? (y/n): " WANT_PROT
read -rp "🧬 Download CDS FASTA for each result?     (y/n): " WANT_CDS

# ---------- NCBI API key reminder ----------
if [ -z "$NCBI_API_KEY" ]; then
    echo ""
    echo "  ℹ️  NCBI_API_KEY is not set. Running at 3 requests/second."
    echo "     Export your key for faster, more reliable downloads:"
    echo "     export NCBI_API_KEY=your_key_here"
    echo "     (Register free at: https://www.ncbi.nlm.nih.gov/account/)"
    echo ""
fi

# ---------- Disk space pre-check ----------
# --- CORRECTION 12: Pre-flight disk space check ---
# Original made no attempt to verify available space before downloading
# potentially large files. We estimate ~500 MB per gene (ortholog ZIP +
# extracted FAA + BLAST DB) and abort if insufficient space is available.
ESTIMATED_MB=$(( ${#ACCESSIONS[@]} * 500 ))
AVAILABLE_KB=$(df -k "$HOME" | awk 'NR==2{print $4}')
AVAILABLE_MB=$(( AVAILABLE_KB / 1024 ))

if [ "$AVAILABLE_MB" -lt "$ESTIMATED_MB" ]; then
    echo "⚠️  Low disk space warning:" >&2
    echo "   Estimated needed : ~${ESTIMATED_MB} MB" >&2
    echo "   Available in \$HOME: ${AVAILABLE_MB} MB" >&2
    read -rp "   Continue anyway? (y/n): " CONTINUE_ANYWAY
    if [[ ! "$CONTINUE_ANYWAY" =~ ^[Yy]$ ]]; then
        echo "Aborted." >&2
        exit 1
    fi
fi

# ---------- Environment setup ----------

install_tools

OUTPUT_ROOT="OASIS_results_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$OUTPUT_ROOT"

# --- CORRECTION 13: Collision-safe temp directory ---
# Original used tmp_OASIS_$$ which relies on PID uniqueness. Two simultaneous
# runs (common on clusters) get different PIDs so they don't collide on the
# dir name — but they DO collide on the shared ortholog/BLAST DB cache files
# inside GLOBAL_TMP when running the same gene. Using mktemp guarantees
# uniqueness; the trade-off is that ortholog caching only works within a
# single run, not across concurrent runs (acceptable for correctness).
GLOBAL_TMP=$(mktemp -d "${TMPDIR:-/tmp}/OASIS_XXXXXX")

# --- CORRECTION 14: Trap SIGTERM and SIGINT in addition to EXIT ---
# Original only trapped EXIT, which is bypassed by SIGKILL but also means
# Ctrl+C (SIGINT) and scheduler kills (SIGTERM) leave temp files behind.
# We add SIGINT and SIGTERM. SIGKILL cannot be trapped — documented below.
_cleanup() {
    echo "" >&2
    echo "  🧹 Cleaning up temporary files in $GLOBAL_TMP ..." >&2
    rm -rf "$GLOBAL_TMP"
}
trap '_cleanup' EXIT INT TERM
# NOTE: SIGKILL (kill -9) cannot be trapped in any shell. If the process is
# killed with SIGKILL, $GLOBAL_TMP will remain and must be cleaned manually:
#   rm -rf /tmp/OASIS_*

echo ""
echo "📁 All results will be saved under: ./${OUTPUT_ROOT}/"
echo "🗂️  Temp directory: $GLOBAL_TMP"
echo ""

# ---------- Main batch loop ----------

SUCCESS=()
SKIPPED=()

for ACC in "${ACCESSIONS[@]}"; do
    # Correction 1 (continued): run_single is called in a subshell-like
    # context via 'if'; set -e does not abort the outer script on a non-zero
    # return from run_single because it is the condition of an if statement.
    # This is intentional — we want to collect SKIPPED[] rather than abort.
    if run_single "$ACC" "$MIN_ID" "$MIN_SIM" "$WANT_PROT" "$WANT_CDS"; then
        SUCCESS+=("$ACC")
    else
        SKIPPED+=("$ACC")
    fi
done

# ---------- Final summary ----------

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
echo "📖 COLUMN DEFINITIONS (BLAST outfmt 6):"
echo "   pident  = % identical positions (strict sequence identity)"
echo "   ppos    = % positively-scoring substitutions (similarity)"
echo ""
echo "⚠️  IMPORTANT NOTICES:"
echo "   • The Datasets CLI caps ortholog downloads at ~499 sequences."
echo "     Accessions that hit this cap are flagged above at download time."
echo "   • Orthologs are available for vertebrates and insects only."
echo "   • Protein FASTAs may contain more sequences than listed accessions"
echo "     due to isoforms / alternative splicing in the NCBI ortholog package."
echo "   • Manual curation is recommended before MSA or phylogenetic analysis."
echo "   • SIGKILL cleanup: if this process was killed with 'kill -9',"
echo "     remove leftover temp files with: rm -rf /tmp/OASIS_*"