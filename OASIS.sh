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

# --- 0. Helper utilities (rootless / sudo-free) ---

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

download_text() {
    local url="$1"
    local result

    for attempt in 1 2 3 4 5; do

        if command -v curl &>/dev/null; then

            result=$(curl \
                --silent \
                --location \
                --fail \
                --retry 3 \
                --retry-delay 2 \
                "$url" 2>/dev/null)

        elif command -v wget &>/dev/null; then

            result=$(wget -q -O - "$url" 2>/dev/null)

        fi

        if [ -n "$result" ]; then
            printf '%s' "$result"
            return 0
        fi

        echo "  ⚠️ Retry $attempt/5..." >&2
        sleep 2
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
#
# Returns one of: protein | nucleotide | gene_id
# protein    → NP_ XP_ WP_  (query with blastp, fetch from protein db)
# nucleotide → NM_ XM_ NG_  (query with blastx, fetch from nuccore db)
# gene_id    → purely numeric strings (treated as NCBI Gene IDs directly)

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
#
# All accession types are normalised to a numeric NCBI Gene ID here.
# protein    → esearch(protein)  → elink(protein→gene)
# nucleotide → esearch(nuccore) → elink(nuccore→gene)
# gene_id    → passed through unchanged

resolve_gene_id() {

    local accession="$1"
    local acc_type
    acc_type=$(detect_type "$accession")

    case "$acc_type" in

    protein)

        echo "  🔎 Resolving protein accession → Gene ID..." >&2

        local puid

        puid=$(download_text \
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=protein&term=${accession}&retmode=json" |
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
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=protein&db=gene&id=${puid}" |
    python3 -c "
import sys
import xml.etree.ElementTree as ET

try:
    root = ET.fromstring(sys.stdin.read())

    accepted = {
        'protein_gene',
        'protein_gene_refseq'
    }

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
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=${accession}&retmode=json" |
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
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=nuccore&db=gene&id=${nuid}" |
    python3 -c "
import sys
import xml.etree.ElementTree as ET

try:
    root = ET.fromstring(sys.stdin.read())

    accepted = {
        'nuccore_gene',
        'nucleotide_gene',
        'nuccore_gene_refseq'
    }

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
#
# Downloads ALL ortholog proteins for a given NCBI Gene ID using the
# Datasets CLI (--ortholog all). The orthologs are only available for vertebrates and insects
# Deduplication of headers is done with awk after concatenation.

fetch_all_orthologs() {
    local gene_id="$1"
    local output_faa="$2"
    local ortho_zip="$GLOBAL_TMP/orthologs_${gene_id}.zip"
    local ortho_dir="$GLOBAL_TMP/orthologs_${gene_id}"

    # Skip if already fetched in this session (shared across queries)
    if [ -s "$output_faa" ]; then
        local cached
        cached=$(grep -c '^>' "$output_faa" 2>/dev/null || echo 0)
        echo "  ♻️  Using cached ortholog FAA ($cached sequences)."
        return 0
    fi

    echo "  🌍 Downloading ortholog protein sequences (gene ID: $gene_id)..."

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

    # Concatenate all protein FAA files found in the package
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

    if (seen[id]++) {
        skip = 1
    } else {
        skip = 0
        print
    }

    next
}

!skip {
    print
}
' "$output_faa" > "${output_faa}.dedup"

mv "${output_faa}.dedup" "$output_faa"

    local count
    count=$(grep -c '^>' "$output_faa" 2>/dev/null || echo 0)
    echo "  ✅ $count ortholog protein sequences ready."
}

# --- 5. Per-query BLAST + filtering ---
#
# Builds a local protein BLAST DB from the ortholog FAA, runs blastp or
# blastx against the query, and filters hits by identity and similarity.

run_blast_filter() {
    local query_fasta="$1" # query sequence file
    local ortho_faa="$2"   # ortholog protein database source
    local blast_prog="$3"  # blastp or blastx
    local min_id="$4"
    local min_sim="$5"
    local db_prefix="$6"  # path prefix for the blast DB files
    local out_list="$7"   # output accession list
    local exclude_id="$8" # accession to exclude (the query itself)

    # Build DB only once per ortholog set (shared across queries of same gene)
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
    "$BLAST_DIR/$blast_prog" \
    -query "$query_fasta" \
    -db "$db_prefix" \
    -max_hsps 1 \
    -max_target_seqs 5000 \
    -outfmt "6 saccver pident ppos" \
    -evalue 1e-5 |
    awk -v id_min="$min_id" -v sim_min="$min_sim" '
        ($2+0) >= (id_min+0) && ($3+0) >= (sim_min+0) {
            print $1
        }
    ' |
    awk -v q="$exclude_id" '
        BEGIN {
            sub(/\.[0-9]+$/, "", q)
        }

        {
            x = $0
            sub(/\.[0-9]+$/, "", x)

            if (x != q)
                print $0
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

    local batch_size=200

    split -l $batch_size "$acclist" "${GLOBAL_TMP}/cds_batch_"

    for batch in "${GLOBAL_TMP}"/cds_batch_*; do

        [ -f "$batch" ] || continue

        local ids
        ids=$(paste -sd, "$batch")

        local tmp_fa="${GLOBAL_TMP}/tmp_cds.fa"

        local success=0

        for retry in {1..5}; do

            download_text \
"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=${ids}&rettype=fasta_cds_na&retmode=text" \
> "$tmp_fa"

            if [[ -s "$tmp_fa" ]]; then
                success=1
                break
            fi

            echo "    ⚠️ Retry $retry/5..."
            sleep 2
        done

        if [[ $success -eq 1 ]]; then

            cat "$tmp_fa" >> "$out_fasta"

            local nseq
            nseq=$(grep -c "^>" "$tmp_fa")

            echo "    → Batch $(basename "$batch"): $nseq CDS sequences"

        else

            echo "    ⚠️ Batch failed."

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
#
# Called once per accession in the batch. All intermediate files go to
# GLOBAL_TMP; final outputs go to OUTPUT_ROOT/<accession>/.

run_single() {
    local acc="$1"
    local min_id="$2"
    local min_sim="$3"
    local want_prot="$4" # y/n
    local want_cds="$5"  # y/n

    echo ""
    echo "╔══════════════════════════════════════════════════════╗"
    echo "║  Processing: $acc"
    echo "╚══════════════════════════════════════════════════════╝"

    # 7a. Detect type
    local acc_type
    acc_type=$(detect_type "$acc")
    if [ "$acc_type" = "unknown" ]; then
        echo "  ⏭️  Skipping '$acc' (unrecognised format)."
        return 1
    fi
    echo "  🔬 Type detected: $acc_type"

    # 7b. Resolve Gene ID
    echo "  🔎 Resolving Gene ID..."
    local gene_id
    gene_id=$(resolve_gene_id "$acc")
    if [ -z "$gene_id" ]; then
        echo "  ❌ Could not resolve Gene ID for '$acc'. Skipping." >&2
        return 1
    fi
    echo "  ✅ Gene ID: $gene_id"

    # 7c. Per-accession output directory
    local out_dir="${OUTPUT_ROOT}/${acc}"
    mkdir -p "$out_dir"

    local final_list="${out_dir}/filtered_accessions_ID${min_id}_SIM${min_sim}_${acc}.txt"

    # 7d. Fetch query sequence
    echo "  📥 Fetching query sequence..."
    local query_fasta="$GLOBAL_TMP/query_${acc}.fasta"
    local blast_prog

    if [ "$acc_type" = "protein" ]; then
        blast_prog="blastp"
        download_text \
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=${acc}&rettype=fasta&retmode=text" \
            >"$query_fasta"
    elif [ "$acc_type" = "nucleotide" ]; then
        blast_prog="blastx"
        download_text \
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${acc}&rettype=fasta&retmode=text" \
            >"$query_fasta"
    else
        # Numeric gene ID: fetch the representative RefSeq protein via esummary
        blast_prog="blastp"
        echo "  ℹ️  Numeric gene ID — fetching representative protein via esummary..."
        local rep_prot
        rep_prot=$(download_text \
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=${acc}&retmode=json" |
            python3 -c "
import sys, json
try:
    data = json.load(sys.stdin)
    uid = list(data['result'].keys() - {'uids'})[0]
    prots = data['result'][uid].get('genomicinfo', [{}])
    # Try to get the accession from protein products
    prod = data['result'][uid].get('locationhist', [{}])
    # fallback: just report the gene_id for manual resolution
    print(data['result'][uid].get('name',''))
except: pass
" 2>/dev/null)
        echo "  ℹ️  Gene symbol: $rep_prot — fetching protein FASTA via gene link..."
        # For gene IDs, fetch the protein via elink gene→protein then efetch
        local prot_uid
        prot_uid=$(download_text \
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=gene&db=protein&id=${acc}&retmode=json" |
            python3 -c "
import sys, json
try:
    data = json.load(sys.stdin)
    for ls in data.get('linksets', []):
        for db in ls.get('linksetdbs', []):
            if db.get('linkname') == 'gene_protein_refseq':
                lk = db.get('links', [])
                if lk: print(lk[0]); raise SystemExit
except: pass
" 2>/dev/null)
        if [ -n "$prot_uid" ]; then
            download_text \
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=${prot_uid}&rettype=fasta&retmode=text" \
                >"$query_fasta"
        fi
    fi

    if [ ! -s "$query_fasta" ] || ! grep -q '^>' "$query_fasta"; then
        echo "  ❌ Could not fetch query sequence for '$acc'. Skipping." >&2
        return 1
    fi
    echo "  ✅ Query sequence fetched."

    # 7e. Fetch orthologs (shared per gene_id across the whole run)
    local ortho_faa="$GLOBAL_TMP/ortho_${gene_id}.faa"
    local blast_db="$GLOBAL_TMP/blastdb_${gene_id}"

    if ! fetch_all_orthologs "$gene_id" "$ortho_faa"; then
        echo "  ❌ Ortholog fetch failed for '$acc'. Skipping." >&2
        return 1
    fi

    # 7f. BLAST + filter
    echo "  ⚙️  Running BLAST alignment and filtering..."
    run_blast_filter \
        "$query_fasta" "$ortho_faa" "$blast_prog" \
        "$min_id" "$min_sim" \
        "$blast_db" "$final_list" "$acc"

    local count
    count=$(wc -l <"$final_list" | tr -d '[:space:]')
    echo "  🎯 $count accessions passed filters (Identity ≥ ${min_id}%, Similarity ≥ ${min_sim}%)."
    echo "  📋 Accession list → $final_list"

    if [ "$count" -eq 0 ]; then
        echo "  ⚠️  No hits passed filters — try lowering thresholds."
        return 0
    fi

    # 7g. Protein FASTA (optional)
    if [[ "$want_prot" =~ ^[Yy]$ ]]; then
        local prot_out="${out_dir}/sequences_PROT_OASIS_${acc}.fasta"
        extract_protein_fasta "$blast_db" "$final_list" "$prot_out"
    fi

    # 7h. CDS FASTA (optional)
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

    ⠀⠀⠉⠓⢦⣄⡀⠀⠉⠙⠲⢼⣧⡉⠙⠲⣤⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
    ⠀⠀⠀⠀⠀⠈⠙⠳⣤⡀⠀⢀⠈⠙⠲⣄⣄⠙⢦⣴⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
    ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⠳⣌⡓⢤⡀⠈⠉⢣⠀⠻⠈⢷⠀⠀⢀⣀⣠⣤⣤⣤⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
    ⣠⣤⡤⠦⠤⠤⠤⠤⠤⠤⢤⣼⡿⣆⠙⢦⠀⠀⢧⠹⡄⠈⡷⠛⢉⣀⣀⡀⠀⠐⠻⠶⠤⢤⣄⡀⠀⠀⠀⠀⠀⠀⠀
    ⠉⠙⠛⠛⠓⠲⠶⠦⢤⣄⡀⠈⢡⡈⠑⢦⡱⡄⠈⣇⠁⢀⠇⠀⣉⡤⠔⢛⣧⡤⠤⠤⠤⠤⠤⠿⢷⣦⣶⣦⣤⣤⠄
    ⠀⠀⠀⠀⠀⠀⠀⠀⢀⣾⣯⣤⣀⡉⠒⠀⣹⠶⠤⣼⣄⣾⠔⠋⠁⣀⡤⠞⢁⣀⣠⠤⠶⠶⠿⣭⣉⠙⢧⡀⠀⠀⠀
    ⠀⠀⠀⠀⣤⣶⡟⠋⠉⣉⣉⣉⠛⠛⠷⣶⡃⢀⡤⠘⡿⠓⢶⣬⠥⠔⠒⠒⠚⠛⠶⢶⣦⠀⠀⠀⠈⠉⠛⠿⠀⠀⠀
    ⠀⠀⠀⢠⡼⠋⠁⣠⣼⣧⣠⣤⠴⠶⠶⠾⢷⣆⣀⣴⠃⢷⣀⡧⢤⣗⡒⠶⠤⠀⠀⠻⢦⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀
    ⠀⠀⡴⠋⢀⡴⠋⠀⠉⠉⠀⠀⠀⠀⠀⢀⡞⠁⠀⢉⡿⠛⣿⢀⢦⡀⠉⢳⣶⠶⠶⢤⣄⡈⠙⠶⣄⠀⠀⠀⠀⠀⠀
    ⠀⣼⠁⡴⠋⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⡎⠀⠀⠀⡼⠀⠀⠛⠻⣆⠳⡀⠘⣏⠁⠀⠀⠀⠉⠙⠓⠮⢿⣦⡀⠀⠀⠀
    ⠀⣇⡞⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡾⠉⠙⠒⢲⡇⠀⠀⠀⠀⠸⡄⢱⠀⢹⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
    ⠸⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠇⢤⠀⠀⢸⠁⠀⠀⠀⠀⠀⢻⠀⠃⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
    ⠀⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢼⡀⠀⠘⠀⢸⠀⠀⠀⠀⠀⠀⢸⡆⢀⡾⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
    ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣾⠉⠉⠉⠉⠙⡇⠀⠀⠀⠀⠀⣸⣧⠞⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
    ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣿⠤⠀⠀⠀⠀⢷⠀⠀⠀⠀⠀⠿⠋⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
    ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢿⠀⠉⠀⠀⠀⠘⡆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
    ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⡀⠀⠀⣀⣀⣠⠼⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
    ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⣿⠋⠉⠁⠀⢀⡀⠻⣆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
    ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠹⡄⠀⠀⠀⠈⠉⠀⠘⢦⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
    ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡀⣀⣀⢹⡀⣰⠛⢧⣠⢄⣀⣬⢧⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
    ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⣿⠻⡌⠻⠿⠃⠀⠈⠁⠈⠁⠸⠋⢹⡷⣶⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
    ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠛⠲⠶⠒⠒⠚⠛⠛⠛⠛⠓⠓⠛⠛⠋⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
    ⠀⠀⠀⠀⠀⠀⠘⠚⠛⠚⠛⠻⠟⠛⠃⠟⠻⠓⠐⠚⠛⠛⠀⠀⠀⠀⠀⠃⠟⠓⠺⠛⠻⠗⠃⠀⠀⠀⠀⠀⠀⠀⠀

          Ortholog Alignment & Similarity Screener
    EvoMol - LAboratório de Evolução Molecular e Sistemas
=============================================================
OASIS_BANNER

# ---------- Input parsing ----------
#
# Accepted call signatures:
#   ./OASIS.sh                             → fully interactive
#   ./OASIS.sh ACC 90 95                   → single accession
#   ./OASIS.sh ACC1 ACC2 ... 90 95         → multiple inline accessions
#   ./OASIS.sh file.txt 90 95             → batch file

ACCESSIONS=()
MIN_ID=""
MIN_SIM=""

if [ "$#" -eq 0 ]; then
    # Fully interactive
    echo ""
    echo "You can enter:"
    echo "  • A single accession:           NP_001416352.1"
    echo "  • Multiple accessions (space):  NP_001416352.1 NM_001429423.1"
    echo "  • A path to a batch file:       /path/to/ids.txt"
    echo ""
    read -rp "🧬 Accession(s) or batch file: " INPUT_RAW

    # Check if it's a file
    if [ -f "$INPUT_RAW" ]; then
        mapfile -t ACCESSIONS < <(grep -v '^\s*#' "$INPUT_RAW" | grep -v '^\s*$' | awk '{print $1}')
    else
        read -ra ACCESSIONS <<<"$INPUT_RAW"
    fi

    read -rp "📊 Minimum Identity (e.g. 90): " MIN_ID
    read -rp "📊 Minimum Similarity (e.g. 95): " MIN_SIM

else
    # CLI mode: last two numeric args are MIN_ID and MIN_SIM
    ARGS=("$@")
    N=${#ARGS[@]}

    # Check if last two look numeric
    if [[ "${ARGS[$((N - 1))]}" =~ ^[0-9]+$ ]] && [[ "${ARGS[$((N - 2))]}" =~ ^[0-9]+$ ]]; then
        MIN_SIM="${ARGS[$((N - 1))]}"
        MIN_ID="${ARGS[$((N - 2))]}"
        ARGS=("${ARGS[@]:0:$((N - 2))}")
    fi

    # Remaining args: single file or list of accessions
    if [ "${#ARGS[@]}" -eq 1 ] && [ -f "${ARGS[0]}" ]; then
        mapfile -t ACCESSIONS < <(grep -v '^\s*#' "${ARGS[0]}" | grep -v '^\s*$' | awk '{print $1}')
        echo "📂 Batch file: ${ARGS[0]} (${#ACCESSIONS[@]} accessions)"
    else
        ACCESSIONS=("${ARGS[@]}")
    fi

    # Prompt for missing thresholds
    [ -z "$MIN_ID" ] && read -rp "📊 Minimum Identity (%):    " MIN_ID
    [ -z "$MIN_SIM" ] && read -rp "📊 Minimum Similarity (%): " MIN_SIM
fi

# Validate we have something to process
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
echo "📊 Identity ≥ ${MIN_ID}%  |  Similarity ≥ ${MIN_SIM}%"
echo ""

# ---------- Output downloads prompt (asked once for the whole batch) ----------

read -rp "📥 Extract protein FASTA for each result? (y/n): " WANT_PROT
read -rp "🧬 Download CDS FASTA for each result?     (y/n): " WANT_CDS

# ---------- Environment setup ----------

install_tools

# Global output root and shared temp directory
OUTPUT_ROOT="OASIS_results_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$OUTPUT_ROOT"

GLOBAL_TMP="tmp_OASIS_$$"
mkdir -p "$GLOBAL_TMP"
trap 'rm -rf "$GLOBAL_TMP"' EXIT

echo ""
echo "📁 All results will be saved under: ./${OUTPUT_ROOT}/"
echo ""

# ---------- Main batch loop ----------

SUCCESS=()
SKIPPED=()

for ACC in "${ACCESSIONS[@]}"; do
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
echo "⚠️  IMPORTANT NOTICE:"
echo "   • The Datasets CLI caps ortholog downloads at ~499 sequences."
echo "     If your filtered list reached 499, additional orthologs may exist."
echo "   • Protein FASTAs may contain more sequences than listed accessions"
echo "     due to isoforms / alternative splicing in the NCBI ortholog package."
echo "   • Manual curation is recommended before MSA or phylogenetic analysis."
