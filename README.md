# OASIS 🌴 (Ortholog Alignment & Similarity Screener)

OASIS is a robust, interactive command-line pipeline designed for bioinformatics researchers to effortlessly fetch, align, and strictly filter orthologous sequences from NCBI.

Unlike web-based BLAST searches that query entire databases (introducing noise from paralogs and synthetic sequences), OASIS uses NCBI's evolutionary curation to download true orthologs and applies rigorous dual-filtering thresholds (Identity and Similarity/Positives) to generate highly reliable datasets for multiple sequence alignments (MSA) and phylogenetic analyses.

---

## ✨ Key Features

- **Strict Dual-Filtering:** Filters orthologs mathematically by both `% Identity` and `% Similarity` simultaneously.
- **Molecule Flexibility:** Automatically handles both Protein (`NP_`, `XP_`) and Nucleotide (`NM_`, `XM_`) queries, adjusting the internal engine (`blastp` or `blastx`) accordingly.
- **Sudo-Free (Rootless) Architecture:** Perfect for university lab environments. It relies on native Linux tools (`curl`/`wget`, `unzip`/`python3`) and installs required bioinformatics tools locally without needing administrator privileges.
- **Zero-Config Dependencies:** Automatically downloads and configures **NCBI BLAST+ (v2.13.0)** and **NCBI Datasets CLI** in the user's `$HOME` directory.
- **Organized Output:** All output files are automatically saved inside a dedicated folder named after the query Accession ID.
- **Batch Extraction:** Instantly generates ready-to-use multifasta files for both Proteins and Coding DNA Sequences (CDS).
- **Clean Execution:** Uses isolated, process-specific temporary directories that self-destruct upon completion or cancellation, leaving your workspace completely clean.

---

## 🚀 Getting Started

### Prerequisites

OASIS is designed to run on Linux/Unix-based systems. It requires an active internet connection to query NCBI databases.

### Installation

Clone this repository and make the script executable:

```bash
git clone https://github.com/RodrigoOrvate/OASIS.git
cd OASIS
chmod +x OASIS.sh
```

### Running the Pipeline

Simply execute the script. It will guide you through an interactive menu:

```bash
./OASIS.sh
```

### Running with Docker

You can also build a lightweight multi-stage image with the NCBI tools prepared during the build stage. Run the build command from inside the `OASIS` project folder, where the `Dockerfile` and `OASIS.sh` files are located:

```bash
cd OASIS
docker build -t oasis .
```

Run it interactively and mount the current directory so the generated files are written back to your workspace:

```bash
docker run --rm -it -v "$PWD:/data" oasis
```

The container uses `/data` as the working directory, so the output files listed below will appear in the mounted folder. You can stop the interactive pipeline at any time with `Ctrl+C`.

After building, Docker may keep the base images used by the multi-stage build (`alpine:3.20` and `debian:bookworm-slim`) as local cache. Keeping them is fine and makes future builds faster, but you can optionally remove them to save disk space:

```bash
docker rmi alpine:3.20 debian:bookworm-slim
```

### Running with Singularity/Apptainer

After the Docker image is published to Docker Hub, you can pull it with Singularity or Apptainer and run it interactively:

```bash
singularity pull oasis.sif docker://rodrigoorvate/oasis:latest
singularity run --bind "$PWD:/data" oasis.sif
```

If your system uses Apptainer, the equivalent commands are:

```bash
apptainer pull oasis.sif docker://rodrigoorvate/oasis:latest
apptainer run --bind "$PWD:/data" oasis.sif
```

---

## 🛠️ Usage Example

When you run OASIS, you will be prompted to provide your parameters:

1. **Accession ID:** Provide a valid NCBI RefSeq ID (e.g., `NP_001416352.1`).
2. **Molecule Type:** Specify if your query is a Protein (`p`) or Nucleotide (`n`).
3. **Thresholds:** Provide the minimum Identity and Similarity percentages (e.g., `90 95`).

**Interactive Prompts:** After the alignment and filtering steps, OASIS will ask if you want to automatically extract the output files:

- *Do you want to extract the protein FASTA file for these sequences? (y/n)*
- *Do you want to download the nucleotide sequences (CDS) for these orthologs? (y/n)*

---

## 📂 Output Files

All output files are saved inside a folder named after the query Accession ID (e.g., `NP_001416352.1/`). Depending on your choices, OASIS will generate up to three files:

```
NP_001416352.1/
├── filtered_accessions_ID{min}_SIM{max}_{ID}.txt
├── sequences_PROT_OASIS_{ID}.fasta
└── sequences_CDS_OASIS_{ID}.fasta
```

1. `filtered_accessions_ID{min}_SIM{max}_{ID}.txt`: A text file containing the accession numbers of the orthologs that passed the dual-filtering step.
2. `sequences_PROT_OASIS_{ID}.fasta`: A multifasta file containing the amino acid sequences for all filtered orthologs.
3. `sequences_CDS_OASIS_{ID}.fasta`: A multifasta file containing the nucleotide Coding DNA Sequences (CDS) for all filtered orthologs.

---

## ⚠️ Known Limitations

> **Please read carefully before using the output files in downstream analyses (MSA, phylogenetics, etc.).**

### 1. Ortholog retrieval cap (~499 sequences)

The `ncbi-datasets-cli` tool currently returns a maximum of approximately **499 sequences** per ortholog download request. If your filtered accession list reaches exactly 499, this likely indicates that additional orthologs exist in the NCBI database but were silently truncated. The total number of available orthologs for your gene of interest may be higher than what OASIS retrieves.

### 2. Isoform redundancy in the Protein FASTA

The ortholog package downloaded from NCBI (`protein.faa`) contains **all annotated protein isoforms** for each orthologous gene, including products of alternative splicing. As a result, the protein FASTA generated by OASIS may contain **more sequences than accessions listed** in the filtered `.txt` file. For example, a filtered list of 499 accessions may yield 600+ protein sequences in the FASTA, because a single gene locus can have multiple associated isoforms in the database.

### 3. Isoforms and alternative splicing variants may cause issues in alignment tools

Some accessions in the protein FASTA may correspond to **alternatively spliced variants** (e.g., shorter isoforms, retained intron products) that differ significantly from the canonical sequence. When loaded into alignment tools such as MEGA, these sequences may:

- **Not appear or not align correctly**, due to length discrepancies or low similarity to the canonical isoform.
- **Introduce noise** into the multiple sequence alignment (MSA), potentially affecting phylogenetic tree topology.

**Recommendation:** After running OASIS, manually inspect the output FASTA files before proceeding to MSA or phylogenetic analysis. Remove or curate sequences that appear to be non-canonical isoforms or that show anomalous alignment behavior.

### 4. Taxonomic scope is limited to Vertebrates and Insects

The `--ortholog all` flag used by OASIS is restricted by NCBI Datasets to **vertebrates and insects only**. Ortholog data for other taxonomic groups (plants, fungi, bacteria, etc.) is not returned by this tool and is therefore outside the scope of OASIS.

### 5. UniProt and AlphaFold IDs are not supported

OASIS relies on `ncbi-datasets-cli`, which only accepts **NCBI RefSeq accessions** (`NP_`, `XP_`, `NM_`, `XM_`). UniProt IDs (e.g., `P09217`) or AlphaFold IDs are not accepted. If you have a UniProt ID, use the [UniProt ID Mapping tool](https://www.uniprot.org/id-mapping) to find the corresponding RefSeq accession before running OASIS.

---

## 🔄 Upcoming NCBI Infrastructure Changes (August 2026)

> **This section is relevant for users and developers who depend on NCBI data retrieval services.**

NCBI has announced significant changes to its data distribution services scheduled for **August 2026**. While the primary announcement concerns the [PMC Article Dataset Distribution Services](https://ncbiinsights.ncbi.nlm.nih.gov/2026/02/12/pmc-article-dataset-distribution-services/), broader changes to NCBI's E-utilities and FTP-based access (including `efetch` via CLI) may also affect automated pipelines.

**What this means for OASIS:**

The current version uses NCBI E-utilities (`efetch`) to download the query sequence. If access to `efetch` via command-line is deprecated or restricted after August 2026, this step may fail silently or return errors.

**OASIS will be updated** to replace `efetch` with a fully `ncbi-datasets-cli`-based retrieval workflow as the August 2026 deadline approaches, ensuring uninterrupted functionality. Watch this repository for updates.

---

## 🗺️ TODO / Roadmap

The following improvements are planned for future versions of OASIS:

- [ ] **Replace `efetch` with `ncbi-datasets-cli`** for query sequence retrieval, to ensure compatibility with NCBI infrastructure changes in August 2026.
- [ ] **Isoform filtering strategy:** Investigate and implement a method to automatically identify and flag (or optionally remove) alternatively spliced isoforms and redundant isoform accessions from the output FASTA files, to reduce noise in downstream MSA and phylogenetic analyses. This is a non-trivial problem given the complexity of splicing annotation in RefSeq.
- [ ] **Lift the 499-sequence cap:** Explore pagination or alternative API strategies to retrieve the full ortholog set when the NCBI Datasets CLI limit is reached.
- [ ] **Taxonomic filter option:** Allow users to specify a taxonomic scope (e.g., only Mammalia, only Actinopterygii) to reduce dataset size and improve relevance.
- [ ] **Accession count mismatch warning:** Add a runtime warning when the number of sequences in the protein FASTA differs from the number of accessions in the filtered list, alerting users to the presence of multiple isoforms per locus.

---

## 🌐 Web Version (No Installation)

Don't have access to a Linux terminal? You can run the entire OASIS pipeline directly in your browser using Google Colab. No installation, configuration, or powerful hardware required.

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/RodrigoOrvate/OASIS/blob/main/OASIS_Colab.ipynb)

---

## ⚖️ License

This project is licensed under the **MIT License**. You are free to use, modify, and distribute this software for academic or commercial purposes, provided that proper credit is given to the original author.

---

## 🏛️ Acknowledgments & Disclaimer

- **NCBI Data & Tools:** OASIS is an independent, open-source wrapper script. It heavily relies on the [NCBI Datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/) and [BLAST+ executables](https://blast.ncbi.nlm.nih.gov/Blast.cgi). This project is **not** officially affiliated with, maintained, or endorsed by the National Center for Biotechnology Information (NCBI) or the National Institutes of Health (NIH).
- **Academic Context:** This tool was developed to support computational biology and bioinformatics research initiatives (PIBIC) at the Federal University of Rio Grande do Norte (UFRN).

---

## 🔬 Developed For

Developed to streamline rigorous ortholog retrieval for phylogenetic and evolutionary conservation analyses.
