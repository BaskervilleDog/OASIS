#!/bin/bash

# =============================================================================
#  OASIS — Ortholog Alignment & Similarity Screener
#  EvoMol - LAboratório de Evolução Molecular e Sistemas
#
#  PURPOSE:
#    Given one or more NCBI accessions, OASIS downloads all known orthologs
#    for the corresponding gene, runs a local BLAST alignment, and filters
#    the results by percent identity (pident) and percent positive-scoring
#    substitutions (ppos). Optional outputs are protein and/or CDS FASTA
#    files for the sequences that pass the thresholds.
#
#  SUPPORTED INPUT FORMATS:
#    NP_ / XP_ / WP_   — RefSeq protein accessions   → blastp
#    NM_ / XM_ / NG_   — RefSeq nucleotide accessions → blastx
#    numeric            — NCBI Gene IDs               → blastp (representative protein)
#
#  USAGE:
#    ./OASIS.sh                              (fully interactive)
#    ./OASIS.sh NP_001416352.1 90 95        (single accession + thresholds)
#    ./OASIS.sh accessions.txt  90 95       (batch file, one accession per line)
#    ./OASIS.sh NP_123.1 NM_456.1 90 95    (multiple accessions inline)
#
#  OUTPUT STRUCTURE (one sub-directory per accession):
#    OASIS_results_YYYYMMDD_HHMMSS/
#    └── <accession>/
#        ├── filtered_accessions_ID<n>_SIM<n>_<acc>.txt   accession list
#        ├── sequences_PROT_OASIS_<acc>.fasta              protein FASTA (optional)
#        └── sequences_CDS_OASIS_<acc>.fasta               CDS FASTA (optional)
#
#  DEPENDENCIES (auto-installed to $HOME on first run):
#    - NCBI Datasets CLI  (ftp.ncbi.nlm.nih.gov)
#    - BLAST+ 2.13.0      (ftp.ncbi.nlm.nih.gov)
#    - curl or wget, unzip or python3, python3 (always required)
#
#  ENVIRONMENT:
#    NCBI_API_KEY   Optional. Raises NCBI rate limit from 3 to 10 req/s.
#                   Register free at https://www.ncbi.nlm.nih.gov/account/
#
#  REFACTORING NOTES (vs original):
#    - resolve_gene_id split into _resolve_protein_to_gene / _resolve_nuccore_to_gene
#    - run_single split into _fetch_query_sequence / _run_outputs
#    - CLI/interactive input split into parse_args / prompt_interactive
#    - Lookup tables (associative arrays) replace repetitive if/elif chains
#    - Each function has cyclomatic complexity (CC) ≤ 6
# =============================================================================

# ---------------------------------------------------------------------------
# Strict mode.
#   -e  exit immediately on any non-zero return (unless in an if condition)
#   -u  treat unset variables as errors
#   -o pipefail  a pipeline fails if any command in it fails, not just the last
# ---------------------------------------------------------------------------
set -euo pipefail

# =============================================================================
# SECTION 0 — NCBI API KEY
# =============================================================================

# Read from the environment if set; default to empty string (no key).
# With a key: 10 req/s.  Without: 3 req/s with possible temporary bans.
NCBI_API_KEY="${NCBI_API_KEY:-}"

# Appends "&api_key=<KEY>" to a URL fragment when a key is configured.
# Called inline inside URL strings: "...${base_url}$(_api_key_param)"
_api_key_param() {
    [ -n "$NCBI_API_KEY" ] && printf '&api_key=%s' "$NCBI_API_KEY" || true
}

# =============================================================================
# SECTION 1 — LOW-LEVEL DOWNLOAD HELPERS
# =============================================================================

# ---------------------------------------------------------------------------
# download_file URL OUTPUT_PATH
#   Downloads a binary or large file directly to disk.
#   Prefers curl; falls back to wget. Hard-exits if neither is available.
#   Used for tool installation (datasets binary, BLAST+ tarball).
# ---------------------------------------------------------------------------
download_file() {
    local url="$1" output="$2"
    if command -v curl &>/dev/null; then
        # --globoff (-g): disables curl's bracket-range URL expansion so that
        # NCBI query qualifiers like [accn] and [gene] reach the server verbatim
        # instead of being interpreted as character ranges (curl error code 3).
        curl -s -L --globoff -o "$output" "$url"
    elif command -v wget &>/dev/null; then
        wget -q -O "$output" "$url"
    else
        echo "❌ Critical Error: Neither 'curl' nor 'wget' found." >&2
        exit 1
    fi
}

# ---------------------------------------------------------------------------
# download_text URL
#   Fetches a URL and prints the response body to stdout.
#   Implements exponential back-off (2 → 4 → 8 → 16 → 32 s) over 5 attempts
#   to handle transient NCBI 429 / 5xx responses.
#   Returns non-zero if all attempts fail.
#
#   CC = 4  (for-loop + curl/wget branch + empty-result check)
# ---------------------------------------------------------------------------
download_text() {
    local url="$1"
    local result wait_secs=2

    for attempt in 1 2 3 4 5; do
        if command -v curl &>/dev/null; then
            # --fail: treat HTTP 4xx/5xx as errors (exit non-zero) so the
            # || true catches them and triggers a retry rather than returning
            # an HTML error page as if it were valid content.
            # --globoff: same bracket-escaping reason as download_file above.
            result=$(curl --silent --location --fail --globoff "$url" 2>/dev/null) || true
        else
            result=$(wget -q -O - "$url" 2>/dev/null) || true
        fi

        # Non-empty response = success; print and return immediately
        if [ -n "$result" ]; then
            printf '%s' "$result"
            return 0
        fi

        echo "  ⚠️  Retry $attempt/5 (waiting ${wait_secs}s)..." >&2
        sleep "$wait_secs"
        wait_secs=$(( wait_secs * 2 ))   # double the wait on every failure
    done

    return 1   # all five attempts exhausted
}

# ---------------------------------------------------------------------------
# extract_zip ZIP_FILE DEST_DIR
#   Extracts a ZIP archive to a destination directory.
#   Prefers the system 'unzip' command; falls back to Python's zipfile module
#   for environments where unzip is not installed.
# ---------------------------------------------------------------------------
extract_zip() {
    local zip_file="$1" dest_dir="$2"
    if command -v unzip &>/dev/null; then
        # -q quiet, -o overwrite without prompting
        unzip -q -o "$zip_file" -d "$dest_dir"
    elif command -v python3 &>/dev/null; then
        python3 -c "import zipfile; zipfile.ZipFile('$zip_file','r').extractall('$dest_dir')"
    else
        echo "❌ Critical Error: Neither 'unzip' nor 'python3' found." >&2
        exit 1
    fi
}

# =============================================================================
# SECTION 2 — TOOL INSTALLATION
# =============================================================================

# Paths for the two external binaries OASIS manages locally.
# Both are installed into $HOME so no root/sudo access is required.
BLAST_DIR="$HOME/ncbi-blast-2.13.0+/bin"
DATASETS_PATH="$HOME/datasets"

# Prepend tool directories to PATH so subsequent calls work without full paths.
export PATH="$BLAST_DIR:$HOME:$PATH"

# ---------------------------------------------------------------------------
# install_tools
#   Installs NCBI Datasets CLI and BLAST+ 2.13.0 if not already present.
#
#   Validation step: after downloading each binary/archive, the file magic
#   (ELF header or gzip header) is checked before use.  This catches the
#   common failure mode where a blocked network returns an HTML error page
#   that curl writes to disk silently (exit code 0) — without the check,
#   that HTML file gets chmod +x'd and later fails with a cryptic exec error.
# ---------------------------------------------------------------------------
install_tools() {
    # Check both that the file is executable AND that it actually runs.
    # A previously-corrupted download (e.g. HTML error page saved as binary)
    # satisfies -x but fails --version, triggering a clean re-download.
    if [ ! -x "$DATASETS_PATH" ] || ! "$DATASETS_PATH" --version &>/dev/null; then
        echo "📦 Installing NCBI Datasets CLI..."
        rm -f "$DATASETS_PATH"   # remove any corrupt previous attempt
        download_file \
            'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets' \
            "$DATASETS_PATH"

        # Validate binary: look for ELF magic bytes (Linux executable header).
        # If the download was blocked, the file will contain HTML instead.
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

        # Validate archive: gzip/tar magic bytes expected.
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

# =============================================================================
# SECTION 3 — ACCESSION-TYPE DETECTION
# =============================================================================

# ---------------------------------------------------------------------------
# detect_type ACCESSION
#   Determines the molecular type of an NCBI accession from its prefix.
#   Prints one of: "protein", "nucleotide", "gene_id", or "unknown".
#
#   Prefix mapping:
#     NP_ XP_ WP_  → protein    (RefSeq protein records)
#     NM_ XM_ NG_  → nucleotide (RefSeq mRNA / genomic records)
#     [0-9]*       → gene_id    (numeric NCBI Gene identifiers)
#
#   CC = 3  (case with 4 arms, one is the error/unknown path)
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

# =============================================================================
# SECTION 4 — GENE ID RESOLUTION
# =============================================================================
#
# All NCBI ortholog downloads require a numeric Gene ID. These helpers convert
# an accession of any supported type into that ID via Entrez E-utilities:
#
#   Protein path:    esearch (protein db) → puid → elink (protein→gene)
#   Nucleotide path: esearch (nuccore db) → nuid → elink (nuccore→gene)
#   Gene ID path:    already numeric, passed through unchanged
#
# The two type-specific resolvers are private (prefixed _) and called only
# through the public resolve_gene_id dispatcher to keep each function small.

# ---------------------------------------------------------------------------
# Shared Python snippet: reads JSON from stdin and prints the first UID in
# esearchresult.idlist. Stored in a variable to avoid duplication between
# _resolve_protein_to_gene and _resolve_nuccore_to_gene.
# ---------------------------------------------------------------------------
_PY_FIRST_ID="
import sys, json
try:
    data = json.load(sys.stdin)
    ids = data['esearchresult']['idlist']
    if ids: print(ids[0])
except Exception: pass
"

# ---------------------------------------------------------------------------
# _py_elink_gene ACCEPTED_LINKNAMES
#   Reads an elink XML response from stdin and prints the first Gene ID found
#   under any of the accepted LinkName values.
#
#   Args:
#     $1  space-separated list of LinkName strings to accept, e.g.:
#         "protein_gene protein_gene_refseq"
#
#   The accepted-names filter avoids picking up unrelated link types
#   (e.g. protein_gene_abstract) that can appear in the same elink response.
# ---------------------------------------------------------------------------
_py_elink_gene() {
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
                    raise SystemExit   # stop after the first valid Gene ID
except Exception: pass
"
}

# ---------------------------------------------------------------------------
# _resolve_protein_to_gene ACCESSION
#   Converts a protein accession (NP_/XP_/WP_) to an NCBI Gene ID.
#   Step 1: esearch protein db by accession → internal protein UID (puid)
#   Step 2: elink protein→gene using the puid → Gene ID
#
#   CC = 3  (two if-guards for empty results + one return)
# ---------------------------------------------------------------------------
_resolve_protein_to_gene() {
    local accession="$1"
    echo "  🔎 Resolving protein accession → Gene ID..." >&2

    # Step 1: translate the versioned accession string to an internal NCBI UID
    local puid
    puid=$(download_text \
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=protein&term=${accession}[accn]&retmode=json$(_api_key_param)" |
        python3 -c "$_PY_FIRST_ID")

    if [ -z "$puid" ]; then
        echo "  ❌ Failed to resolve protein UID." >&2
        return 1
    fi

    # Step 2: follow the protein→gene link to get the Gene ID
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

# ---------------------------------------------------------------------------
# _resolve_nuccore_to_gene ACCESSION
#   Converts a nucleotide accession (NM_/XM_/NG_) to an NCBI Gene ID.
#   Step 1: esearch nuccore db by accession → internal nuccore UID (nuid)
#   Step 2: elink nuccore→gene using the nuid → Gene ID
#
#   Mirrors _resolve_protein_to_gene but queries the nuccore database and
#   uses the nuccore-specific link names accepted by _py_elink_gene.
#
#   CC = 3
# ---------------------------------------------------------------------------
_resolve_nuccore_to_gene() {
    local accession="$1"
    echo "  🔎 Resolving nucleotide accession → Gene ID..." >&2

    # Step 1: get internal nuccore UID for this accession
    local nuid
    nuid=$(download_text \
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=${accession}[accn]&retmode=json$(_api_key_param)" |
        python3 -c "$_PY_FIRST_ID")

    if [ -z "$nuid" ]; then
        echo "  ❌ Failed to resolve nuccore UID." >&2
        return 1
    fi

    # Step 2: follow the nuccore→gene link
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

# ---------------------------------------------------------------------------
# resolve_gene_id ACCESSION
#   Public dispatcher: routes any supported accession to the correct resolver.
#   Numeric Gene IDs are already in the target format and pass through as-is.
#
#   CC = 4  (case with 3 active branches + fallthrough error)
# ---------------------------------------------------------------------------
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

# =============================================================================
# SECTION 5 — ORTHOLOG FETCHER
# =============================================================================

# ---------------------------------------------------------------------------
# fetch_all_orthologs GENE_ID OUTPUT_FAA
#   Downloads the complete set of orthologous protein sequences for a gene
#   from the NCBI Datasets API, extracts them, and deduplicates by accession.
#
#   Args:
#     $1  GENE_ID      numeric NCBI Gene identifier
#     $2  OUTPUT_FAA   path where the merged, deduplicated FASTA will be written
#
#   Caching: if OUTPUT_FAA already exists and is non-empty, the download is
#   skipped (the same gene may appear in multiple accessions in a batch run).
#
#   Scope warning: NCBI Datasets ortholog downloads cover vertebrates and
#   insects only. Plant, fungal, and other gene IDs will yield 0 orthologs.
#
#   Cap warning: the Datasets API silently truncates results at ~499 sequences.
#   If that limit is detected, a prominent warning is printed so the user can
#   switch to the NCBI web portal or split the query by taxon group.
#
#   CC = 5  (cache check + zip check + faa check + dedup + cap check)
# ---------------------------------------------------------------------------
fetch_all_orthologs() {
    local gene_id="$1" output_faa="$2"
    local ortho_zip="$GLOBAL_TMP/orthologs_${gene_id}.zip"
    local ortho_dir="$GLOBAL_TMP/orthologs_${gene_id}"

    # Cache hit: re-use FAA from an earlier run within the same OASIS session
    if [ -s "$output_faa" ]; then
        local cached
        cached=$(grep -c '^>' "$output_faa" 2>/dev/null || echo 0)
        echo "  ♻️  Using cached ortholog FAA ($cached sequences)."
        return 0
    fi

    echo "  ⚠️  Note: ortholog downloads cover vertebrates and insects only." >&2
    echo "  🌍 Downloading ortholog protein sequences (gene ID: $gene_id)..." >&2

    # stdout suppressed (progress bars); stderr kept so network errors surface.
    # The --ortholog all flag requests orthologs across all available taxa.
    "$DATASETS_PATH" download gene gene-id "$gene_id" \
        --ortholog all --include protein --filename "$ortho_zip" >/dev/null

    if [ ! -s "$ortho_zip" ]; then
        echo "  ❌ Ortholog download failed for gene ID $gene_id." >&2
        return 1
    fi

    mkdir -p "$ortho_dir"
    extract_zip "$ortho_zip" "$ortho_dir"

    # Merge all FAA files from the extracted archive into one file.
    # The archive may contain one file per taxon group.
    find "$ortho_dir" -type f \( -name "*.faa" -o -name "*.protein.faa" \) \
        -exec cat {} + >"$output_faa"

    if [ ! -s "$output_faa" ]; then
        echo "  ❌ No protein sequences found after extraction." >&2
        return 1
    fi

    # Deduplicate: some orthologs appear in multiple taxon FAA files.
    # Strategy: track the first field of each FASTA header (the accession);
    # keep the first occurrence and skip subsequent ones.
    awk '/^>/{
        header=$0; split(header,a," "); id=a[1]
        if (seen[id]++) { skip=1 } else { skip=0; print }
        next
    } !skip { print }' "$output_faa" > "${output_faa}.dedup"
    mv "${output_faa}.dedup" "$output_faa"

    local count
    count=$(grep -c '^>' "$output_faa" 2>/dev/null || echo 0)
    echo "  ✅ $count ortholog protein sequences ready."

    # Warn prominently when the NCBI Datasets hard cap (~499) is hit.
    # Exact count of 499 strongly suggests truncation; results will be incomplete.
    if [ "$count" -ge 499 ]; then
        echo "" >&2
        echo "  ⚠️  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━" >&2
        echo "  ⚠️  DATASET CAP REACHED: $count sequences (NCBI cap ~499)." >&2
        echo "  ⚠️  Results are INCOMPLETE. Split by taxon or use the web portal." >&2
        echo "  ⚠️  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━" >&2
        echo "" >&2
    fi
}

# =============================================================================
# SECTION 6 — BLAST FILTER
# =============================================================================

# ---------------------------------------------------------------------------
# run_blast_filter QUERY_FASTA ORTHO_FAA BLAST_PROG MIN_ID MIN_SIM \
#                  DB_PREFIX OUT_LIST EXCLUDE_ID
#   Builds a local BLAST protein database from the ortholog FAA (or reuses
#   one built in an earlier step for the same gene), runs an alignment, and
#   writes the accessions of hits that pass both identity thresholds to a
#   plain-text list file.
#
#   Args:
#     $1 QUERY_FASTA   FASTA of the input query sequence
#     $2 ORTHO_FAA     FASTA of all ortholog sequences (the search database)
#     $3 BLAST_PROG    blastp (protein query) or blastx (nucleotide query)
#     $4 MIN_ID        minimum pident — % identical positions (strict)
#     $5 MIN_SIM       minimum ppos  — % positively-scoring substitutions
#     $6 DB_PREFIX     path prefix for makeblastdb output files
#     $7 OUT_LIST      path for the filtered accession list (one per line)
#     $8 EXCLUDE_ID    query accession to remove from results (self-match)
#
#   BLAST output format (-outfmt "6 saccver pident ppos"):
#     col 1  saccver — subject accession with version
#     col 2  pident  — % identical positions
#     col 3  ppos    — % positively-scoring positions (similarity)
#
#   The two awk filters in the pipeline:
#     Filter 1: numeric threshold check on pident (col 2) and ppos (col 3)
#     Filter 2: strip the query's own accession from results (version-agnostic,
#               so NP_001.1 and NP_001.2 are both excluded when query is NP_001.1)
#
#   CC = 3  (db-exists check + two awk conditions counted as one decision each)
# ---------------------------------------------------------------------------
run_blast_filter() {
    local query_fasta="$1" ortho_faa="$2" blast_prog="$3"
    local min_id="$4" min_sim="$5" db_prefix="$6" out_list="$7" exclude_id="$8"

    # Build the database only once per gene; subsequent accessions sharing the
    # same Gene ID (e.g. isoforms) will reuse it, saving significant time.
    if [ ! -f "${db_prefix}.pin" ] && [ ! -f "${db_prefix}.pdb" ]; then
        echo "  🔨 Building local BLAST database..."
        # -parse_seqids enables per-accession retrieval via blastdbcmd later
        "$BLAST_DIR/makeblastdb" -in "$ortho_faa" -dbtype prot \
            -out "$db_prefix" -parse_seqids -logfile /dev/null
    else
        echo "  ♻️  Reusing existing BLAST database."
    fi

    echo "  🔬 Running $blast_prog alignment..."
    echo "      Filters: pident ≥ ${min_id}%  |  ppos ≥ ${min_sim}%"

    "$BLAST_DIR/$blast_prog" \
        -query "$query_fasta" -db "$db_prefix" \
        -max_hsps 1 \           # keep only the best HSP per subject sequence
        -max_target_seqs 5000 \ # collect up to 5000 candidate hits before filtering
        -outfmt "6 saccver pident ppos" \
        -evalue 1e-5 |          # discard very weak hits before threshold checks

    # Filter 1: both pident and ppos must meet their respective minimums.
    # Arithmetic comparison (+0) coerces strings to numbers safely in awk.
    awk -v id_min="$min_id" -v sim_min="$min_sim" \
        '($2+0)>=(id_min+0) && ($3+0)>=(sim_min+0) {print $1}' |

    # Filter 2: remove the query accession itself (self-match).
    # Version suffix is stripped from both the query and each hit before
    # comparison so NP_001.1 matches NP_001.2 and is correctly excluded.
    awk -v q="$exclude_id" '
        BEGIN { sub(/\.[0-9]+$/,"",q) }
        { x=$0; sub(/\.[0-9]+$/,"",x); if(x!=q) print $0 }
    ' | sort -u > "$out_list"   # deduplicate and sort for reproducibility
}

# =============================================================================
# SECTION 7 — OUTPUT EXTRACTORS
# =============================================================================

# ---------------------------------------------------------------------------
# extract_protein_fasta DB_PREFIX ACCLIST OUT_FASTA
#   Retrieves the protein sequences for all accessions in ACCLIST from the
#   local BLAST database and writes them to OUT_FASTA.
#   Uses blastdbcmd with -entry_batch for efficient bulk retrieval.
#
#   CC = 2  (one if-check on output file size)
# ---------------------------------------------------------------------------
extract_protein_fasta() {
    local db_prefix="$1" acclist="$2" out_fasta="$3"
    echo "  🚀 Extracting protein sequences from local database..."

    # -entry_batch accepts a file with one accession per line
    "$BLAST_DIR/blastdbcmd" -db "$db_prefix" -entry_batch "$acclist" \
        -out "$out_fasta" 2>/dev/null   # suppress "Entry not found" warnings

    if [ -s "$out_fasta" ]; then
        local n; n=$(grep -c '^>' "$out_fasta")
        echo "  ✅ Protein FASTA: $n sequences → $(basename "$out_fasta")"
    else
        echo "  ❌ Protein extraction returned empty output." >&2
    fi
}

# ---------------------------------------------------------------------------
# extract_cds_fasta ACCLIST OUT_FASTA
#   Downloads the CDS (coding DNA sequence) for every accession in ACCLIST
#   from NCBI Entrez in batches of 100, appending results to OUT_FASTA.
#
#   Batching is necessary because Entrez has a URL-length limit and also
#   performs better with moderate batch sizes than with thousands of IDs
#   in a single request.
#
#   CC = 4  (for-loop + file-guard + per-batch retry delegated to _fetch_cds_batch)
# ---------------------------------------------------------------------------
extract_cds_fasta() {
    local acclist="$1" out_fasta="$2"
    echo "  🚀 Downloading CDS sequences..."
    > "$out_fasta"   # truncate/create the output file

    # Create a private temp directory for batch files to avoid collisions
    # with other accessions processed in the same OASIS run.
    local cds_batch_dir
    cds_batch_dir=$(mktemp -d "${GLOBAL_TMP}/cds_batches_XXXXXX")

    # Split the accession list into files of 100 lines each (batch_aa, batch_ab, …)
    split -l 100 "$acclist" "${cds_batch_dir}/batch_"

    for batch in "${cds_batch_dir}"/batch_*; do
        [ -f "$batch" ] || continue   # skip if glob matched nothing
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

# ---------------------------------------------------------------------------
# _fetch_cds_batch BATCH_FILE OUT_FASTA
#   Downloads CDS sequences for a single batch of accessions (≤100) from
#   NCBI Entrez using the efetch rettype=fasta_cds_na endpoint.
#   Appends successful results to OUT_FASTA.
#   Implements exponential back-off (same pattern as download_text).
#
#   CC = 3  (retry for-loop + non-empty-file check + implicit return on success)
# ---------------------------------------------------------------------------
_fetch_cds_batch() {
    local batch="$1" out_fasta="$2"
    # Entrez efetch accepts a comma-separated list of IDs
    local ids; ids=$(paste -sd, "$batch")
    local tmp_fa="${batch}.fa"
    local wait_secs=2

    for retry in 1 2 3 4 5; do
        # rettype=fasta_cds_na returns nucleotide CDS sequences in FASTA format
        # keyed on the protein accession (same IDs as in our accession list)
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
    # Non-fatal: the outer loop continues with the next batch
}

# =============================================================================
# SECTION 8 — QUERY SEQUENCE FETCHER
# =============================================================================

# ---------------------------------------------------------------------------
# _fetch_query_sequence ACCESSION ACC_TYPE QUERY_FASTA
#   Downloads the FASTA sequence for the query accession and writes it to
#   QUERY_FASTA. Prints the name of the appropriate BLAST program to stdout
#   so the caller can capture it: blast_prog=$(_fetch_query_sequence ...)
#
#   The mapping from accession type to (BLAST program, Entrez database) is
#   encoded in associative arrays (lightweight lookup tables) rather than
#   if/elif chains, keeping the CC low and making future additions trivial.
#
#   For numeric Gene IDs there is no direct sequence record; the function
#   resolves the representative RefSeq protein via gene→protein_refseq elink
#   and fetches that protein's FASTA instead.
#
#   CC = 4  (gene_id branch + prot_uid empty-check + else fetch + return)
# ---------------------------------------------------------------------------
_fetch_query_sequence() {
    local acc="$1" acc_type="$2" query_fasta="$3"

    # Lookup tables: map accession type to BLAST program and Entrez database.
    # Extending to a new type only requires adding one entry to each array.
    local -A BLAST_PROG=( [protein]=blastp [nucleotide]=blastx [gene_id]=blastp )
    local -A FETCH_DB=(   [protein]=protein [nucleotide]=nuccore )

    local blast_prog="${BLAST_PROG[$acc_type]}"

    if [ "$acc_type" = "gene_id" ]; then
        # Numeric Gene IDs have no direct sequence record in Entrez.
        # Use gene→protein_refseq elink to find the representative protein UID,
        # then fetch that protein's FASTA.
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
        # Protein and nucleotide accessions can be fetched directly from Entrez
        local db="${FETCH_DB[$acc_type]}"
        download_text \
"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=${db}&id=${acc}&rettype=fasta&retmode=text$(_api_key_param)" \
            >"$query_fasta"
    fi

    # Return the blast program name to the caller via stdout
    echo "$blast_prog"
}

# =============================================================================
# SECTION 9 — OUTPUT STAGE
# =============================================================================

# ---------------------------------------------------------------------------
# _run_outputs DB_PREFIX ACCESSION OUT_DIR FINAL_LIST WANT_PROT WANT_CDS
#   Conditionally extracts protein and/or CDS FASTA files based on user choice.
#   Extracted from run_single to keep that function's CC low.
#
#   Args:
#     $1 DB_PREFIX    path prefix of the BLAST database (for blastdbcmd)
#     $2 ACCESSION    original query accession (used in output filenames)
#     $3 OUT_DIR      per-accession output directory
#     $4 FINAL_LIST   path to the filtered accession list
#     $5 WANT_PROT    "y"/"Y" to extract protein FASTA
#     $6 WANT_CDS     "y"/"Y" to download CDS FASTA
#
#   CC = 3  (two independent if-checks + implicit return)
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

# =============================================================================
# SECTION 10 — SINGLE-ACCESSION PIPELINE
# =============================================================================

# ---------------------------------------------------------------------------
# run_single ACCESSION MIN_ID MIN_SIM WANT_PROT WANT_CDS
#   Orchestrates the full analysis pipeline for one accession:
#     1. Detect accession type
#     2. Resolve to a numeric Gene ID
#     3. Fetch the query sequence FASTA
#     4. Download all ortholog protein sequences
#     5. Build a local BLAST database and run the alignment
#     6. Apply pident/ppos filters and write the accession list
#     7. Optionally extract protein and/or CDS FASTA outputs
#
#   Each step is delegated to a dedicated function; run_single itself
#   contains only sequencing logic and early-exit error guards.
#   Returns 0 on success (even if 0 hits pass the filter), 1 on any fatal error.
#
#   CC = 5  (unknown-type check + gene_id check + query check + ortho check
#            + zero-hits check)
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

    # Resolve accession to Gene ID (required for the Datasets ortholog download)
    local gene_id
    gene_id=$(resolve_gene_id "$acc") || {
        echo "  ❌ Could not resolve Gene ID for '$acc'. Skipping." >&2
        return 1
    }
    echo "  ✅ Gene ID: $gene_id"

    # Create a dedicated output directory for this accession
    local out_dir="${OUTPUT_ROOT}/${acc}"
    mkdir -p "$out_dir"

    # Fetch the query sequence; blast_prog (blastp or blastx) is returned
    # on stdout so the correct BLAST flavour is used for each accession type
    local query_fasta="$GLOBAL_TMP/query_${acc}.fasta"
    local blast_prog
    blast_prog=$(_fetch_query_sequence "$acc" "$acc_type" "$query_fasta") || return 1

    # Sanity-check: the downloaded FASTA must be non-empty and contain a header
    if [ ! -s "$query_fasta" ] || ! grep -q '^>' "$query_fasta"; then
        echo "  ❌ Could not fetch query sequence for '$acc'. Skipping." >&2
        return 1
    fi
    echo "  ✅ Query sequence fetched."

    # Ortholog FAA and BLAST DB are keyed on gene_id so multiple accessions
    # sharing the same gene reuse the same files within one OASIS session
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
        return 0   # not a failure; results are simply empty
    fi

    _run_outputs "$blast_db" "$acc" "$out_dir" "$final_list" "$want_prot" "$want_cds"
    echo "  ✅ Done → ./${out_dir}/"
}

# =============================================================================
# SECTION 11 — INPUT PARSING
# =============================================================================

# ---------------------------------------------------------------------------
# prompt_interactive
#   Collects accessions and thresholds interactively when OASIS is run with
#   no command-line arguments. Writes to the global variables ACCESSIONS,
#   MIN_ID, and MIN_SIM.
#
#   Accepts either:
#     • a path to a batch file (one accession per line, # comments ignored)
#     • one or more accessions typed directly (space-separated)
#
#   CC = 2  (one if/else for file-vs-inline input)
# ---------------------------------------------------------------------------
prompt_interactive() {
    echo ""
    echo "You can enter:"
    echo "  • A single accession:           NP_001416352.1"
    echo "  • Multiple accessions (space):  NP_001416352.1 NM_001429423.1"
    echo "  • A path to a batch file:       /path/to/ids.txt"
    echo ""
    read -rp "🧬 Accession(s) or batch file: " INPUT_RAW

    if [ -f "$INPUT_RAW" ]; then
        # Batch file: strip comment lines (#) and blank lines; take first field only
        mapfile -t ACCESSIONS < <(grep -v '^\s*#' "$INPUT_RAW" | grep -v '^\s*$' | awk '{print $1}')
    else
        # Inline: split on whitespace into the ACCESSIONS array
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

# ---------------------------------------------------------------------------
# parse_args ACCS_NAMEREF ID_NAMEREF SIM_NAMEREF [ARGS...]
#   Parses command-line arguments and populates the caller's variables via
#   bash namerefs (local -n), avoiding global variable coupling.
#
#   Argument formats accepted:
#     acc1 [acc2 ...] MIN_ID MIN_SIM      inline accessions + thresholds
#     batch_file.txt  MIN_ID MIN_SIM      batch file + thresholds
#     acc1 [acc2 ...]                     accessions only (prompts for thresholds)
#
#   The last two numeric arguments are always interpreted as MIN_ID and MIN_SIM.
#   Any remaining arguments are treated as accessions or a single batch file.
#
#   CC = 4  (threshold-pair check + batch-file check + two missing-threshold prompts)
# ---------------------------------------------------------------------------
parse_args() {
    local -n _accs="$1"   # nameref: caller's accession array
    local -n _id="$2"     # nameref: caller's MIN_ID variable
    local -n _sim="$3"    # nameref: caller's MIN_SIM variable

    local args=("${@:4}")   # all positional args after the three nameref names
    local n=${#args[@]}

    # If the last two arguments are both integers, treat them as thresholds
    # and remove them from the accession list
    if [[ "${args[$((n-1))]}" =~ ^[0-9]+$ ]] && \
       [[ "${args[$((n-2))]}" =~ ^[0-9]+$ ]]; then
        _sim="${args[$((n-1))]}"
        _id="${args[$((n-2))]}"
        args=("${args[@]:0:$((n-2))}")
    fi

    # Determine whether the remaining argument(s) form a batch file or inline list
    if [ "${#args[@]}" -eq 1 ] && [ -f "${args[0]}" ]; then
        mapfile -t _accs < <(grep -v '^\s*#' "${args[0]}" | grep -v '^\s*$' | awk '{print $1}')
        echo "📂 Batch file: ${args[0]} (${#_accs[@]} accessions)"
    else
        _accs=("${args[@]}")
    fi

    # Prompt for any threshold that was not supplied on the command line
    if [ -z "${_id:-}" ]; then
        echo "  pident = % identical positions  |  ppos = % positively-scoring substitutions"
        read -rp "📊 Minimum pident / Identity  (%): " _id
    fi
    if [ -z "${_sim:-}" ]; then
        read -rp "📊 Minimum ppos   / Similarity (%): " _sim
    fi
}

# =============================================================================
# MAIN
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

# ---------------------------------------------------------------------------
# Input collection
# Delegates entirely to prompt_interactive (no args) or parse_args (with args)
# so the main block itself has minimal branching (CC ≈ 4).
# ---------------------------------------------------------------------------
ACCESSIONS=()
MIN_ID=""
MIN_SIM=""

if [ "$#" -eq 0 ]; then
    prompt_interactive
else
    parse_args ACCESSIONS MIN_ID MIN_SIM "$@"
fi

# Validate that we have at least one accession and both thresholds
[ "${#ACCESSIONS[@]}" -eq 0 ] && { echo "❌ No accessions provided. Exiting." >&2; exit 1; }
[ -z "$MIN_ID" ] || [ -z "$MIN_SIM" ] && { echo "❌ Thresholds required." >&2; exit 1; }

echo ""
echo "📋 Accessions to process (${#ACCESSIONS[@]}):"
for a in "${ACCESSIONS[@]}"; do echo "   • $a"; done
echo "📊 pident ≥ ${MIN_ID}%  |  ppos ≥ ${MIN_SIM}%"
echo ""

read -rp "📥 Extract protein FASTA for each result? (y/n): " WANT_PROT
read -rp "🧬 Download CDS FASTA for each result?     (y/n): " WANT_CDS

# ---------------------------------------------------------------------------
# NCBI API key reminder
# ---------------------------------------------------------------------------
if [ -z "$NCBI_API_KEY" ]; then
    echo ""
    echo "  ℹ️  NCBI_API_KEY not set — running at 3 req/s."
    echo "     export NCBI_API_KEY=your_key_here"
    echo "     (Register free at: https://www.ncbi.nlm.nih.gov/account/)"
    echo ""
fi

# ---------------------------------------------------------------------------
# Disk space pre-check
# Estimate ~500 MB per accession (ortholog ZIP + extracted FAA + BLAST DB).
# Warn and offer to abort if insufficient space is detected.
# ---------------------------------------------------------------------------
ESTIMATED_MB=$(( ${#ACCESSIONS[@]} * 500 ))
AVAILABLE_MB=$(( $(df -k "$HOME" | awk 'NR==2{print $4}') / 1024 ))

if [ "$AVAILABLE_MB" -lt "$ESTIMATED_MB" ]; then
    echo "⚠️  Low disk space: need ~${ESTIMATED_MB} MB, have ${AVAILABLE_MB} MB." >&2
    read -rp "   Continue anyway? (y/n): " CONTINUE_ANYWAY
    [[ ! "$CONTINUE_ANYWAY" =~ ^[Yy]$ ]] && { echo "Aborted." >&2; exit 1; }
fi

# ---------------------------------------------------------------------------
# Environment setup
# ---------------------------------------------------------------------------
install_tools

# Timestamped output root keeps successive runs from overwriting each other
OUTPUT_ROOT="OASIS_results_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$OUTPUT_ROOT"

# Collision-safe temp directory; mktemp guarantees uniqueness even across
# concurrent OASIS invocations processing different genes on the same machine
GLOBAL_TMP=$(mktemp -d "${TMPDIR:-/tmp}/OASIS_XXXXXX")

# ---------------------------------------------------------------------------
# Cleanup trap
# Removes the temp directory on normal exit, Ctrl+C (SIGINT), or SIGTERM.
# SIGKILL (kill -9) cannot be trapped in bash; leftover files from a kill -9
# must be removed manually: rm -rf /tmp/OASIS_*
# ---------------------------------------------------------------------------
_cleanup() {
    echo "" >&2
    echo "  🧹 Cleaning up temporary files in $GLOBAL_TMP ..." >&2
    rm -rf "$GLOBAL_TMP"
}
trap '_cleanup' EXIT INT TERM

echo ""
echo "📁 Results: ./${OUTPUT_ROOT}/"
echo "🗂️  Temp:    $GLOBAL_TMP"
echo ""

# ---------------------------------------------------------------------------
# Main batch loop
# Each accession is processed independently; failures are collected in
# SKIPPED[] rather than aborting the entire run. The if-condition around
# run_single is intentional: it prevents set -e from killing the script
# on a per-accession failure.
# ---------------------------------------------------------------------------
SUCCESS=()
SKIPPED=()

for ACC in "${ACCESSIONS[@]}"; do
    if run_single "$ACC" "$MIN_ID" "$MIN_SIM" "$WANT_PROT" "$WANT_CDS"; then
        SUCCESS+=("$ACC")
    else
        SKIPPED+=("$ACC")
    fi
done

# ---------------------------------------------------------------------------
# Final summary
# ---------------------------------------------------------------------------
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
echo "   pident = % identical positions (strict sequence identity)"
echo "   ppos   = % positively-scoring substitutions (similarity)"
echo ""
echo "⚠️  NOTICES:"
echo "   • NCBI caps ortholog downloads at ~499 sequences."
echo "   • Orthologs available for vertebrates and insects only."
echo "   • Manual curation recommended before MSA or phylogenetic analysis."
echo "   • SIGKILL cleanup: rm -rf /tmp/OASIS_*"