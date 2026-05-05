# notebook-log.md
# PLPATH 563 – Phylogenetic Analysis of Molecular Data
**Student:** Aidan Anguyse
**Date Started:** 2026-02-20
**GitHub Repo:** Documents/GitHub/myProject/myProject

---

## Table of Contents
1. [Dataset Description](#dataset-description)
2. [Project Folder Structure](#project-folder-structure)
3. [Software & Tools Used](#software--tools-used)
4. [Step 1: Environment Setup](#step-1-environment-setup)
5. [Step 2: Download Raw Reads](#step-2-download-raw-reads)
6. [Step 3: Compress Raw Reads](#step-3-compress-raw-reads)
7. [Step 4: Quality Control with FastQC](#step-4-quality-control-with-fastqc)
8. [Sample Summary Table](#sample-summary-table)
9. [QC Results & Observations](#qc-results--observations)
10. [Next Steps](#next-steps)

---

## Dataset Description

**Paper:** CRyPTIC Consortium. "A data compendium associating the genomes of 12,289
Mycobacterium tuberculosis isolates with quantitative resistance phenotypes to 13 antibiotics."
*PLOS Biology*, 2022. DOI: 10.1371/journal.pbio.3001721. PMID: 35944069

**BioProject:** PRJEB35849 (European Nucleotide Archive / NCBI SRA)
**Organism:** *Mycobacterium tuberculosis* — the bacterium responsible for tuberculosis (TB)
**Data type:** Illumina paired-end whole-genome sequencing (raw reads, .fastq.gz)
**Full dataset scope:** 12,289 clinical isolates from 23 countries worldwide
**Subset used in this project:** 8 isolates selected to represent phylogenetic and resistance diversity

**Why this dataset:**
This compendium links bacterial genomic variants to quantitative antibiotic resistance
(minimum inhibitory concentrations, MICs) for 13 TB drugs. It spans all 4 major
*M. tuberculosis* lineages (L1–L4) and the full spectrum of resistance — from fully
drug-susceptible to extensively drug-resistant (XDR). This makes it ideal for
phylogenetic analysis of resistance evolution in a globally distributed pathogen.

---

## Project Folder Structure

```
myProject/
├── notebook-log.md               # This file — reproducibility log
├── README.md                     # Project overview
├── .gitignore                    # Tells Git to ignore large raw data files
├── phylogenetics-class/          # Class notes and cheatsheets
│   └── software-cheatsheet.md   # Software descriptions table
├── raw_reads/                    # Raw Illumina .fastq.gz files (NOT pushed to GitHub)
├── qc_results/                   # FastQC .html and .zip output reports
├── data/                         # Processed/downstream data (VCFs, alignments, trees)
└── scripts/                      # Bash scripts for reproducibility
    └── 01_download_and_qc.sh    # Full download + QC pipeline script
```

---

## Software & Tools Used

| Software | Version | Purpose |
|---|---|---|
| SRA Toolkit (`fasterq-dump`) | 3.3.0 | Download raw reads from NCBI SRA |
| FastQC | 0.12.1 | Quality control of raw Illumina reads |
| gzip | system | Compress .fastq files to .fastq.gz |

---

## Step 1: Environment Setup

### 1a. Install SRA Toolkit (Mac ARM64 / Apple Silicon)
The standard `conda install -c bioconda sra-tools` fails on `osx-arm64` due to a
`perl-data-dumper` dependency conflict. Instead, download the binary directly from NCBI:

```bash
# Download the Mac ARM64 binary
curl -o sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-mac-arm64.tar.gz

# Extract
tar -xzf sratoolkit.tar.gz

# Check the folder name that was created
ls
# You should see: sratoolkit.3.3.0-mac-arm64

# Add to PATH
export PATH=$PATH:~/sratoolkit.3.3.0-mac-arm64/bin

# Make PATH permanent across terminal sessions
echo 'export PATH=$PATH:~/sratoolkit.3.3.0-mac-arm64/bin' >> ~/.zshrc
source ~/.zshrc

# Configure SRA toolkit (press 'x' to exit the menu)
vdb-config --interactive

# Verify installation
fasterq-dump --version
# Expected output: fasterq-dump : 3.3.0
```

### 1b. Verify FastQC is installed
FastQC was pre-installed. Confirm it is accessible:
```bash
fastqc --version
```

### 1c. Create project subfolder structure
```bash
mkdir -p ~/Documents/GitHub/myProject/myProject/raw_reads
mkdir -p ~/Documents/GitHub/myProject/myProject/qc_results
mkdir -p ~/Documents/GitHub/myProject/myProject/data
mkdir -p ~/Documents/GitHub/myProject/myProject/scripts
```

---

## Step 2: Download Raw Reads

Navigate to the raw_reads folder first, then download each sample one at a time.
Each download takes approximately 5–15 minutes depending on internet speed.
The `--split-files` flag separates paired-end reads into _1.fastq (forward) and _2.fastq (reverse).

```bash
cd ~/Documents/GitHub/myProject/myProject/raw_reads

# Sample 1 — ERR4796930
fasterq-dump ERR4796930 --split-files --threads 4

# Sample 2 — ERR4796931
fasterq-dump ERR4796931 --split-files --threads 4

# Sample 3 — ERR4796932
fasterq-dump ERR4796932 --split-files --threads 4

# Sample 4 — ERR4796935
fasterq-dump ERR4796935 --split-files --threads 4

# Sample 5 — ERR4796940
fasterq-dump ERR4796940 --split-files --threads 4

# Sample 6 — ERR4796945
fasterq-dump ERR4796945 --split-files --threads 4

# Sample 7 — ERR4796950
fasterq-dump ERR4796950 --split-files --threads 4

# Sample 8 — ERR4796960
fasterq-dump ERR4796960 --split-files --threads 4
```

---

## Step 3: Compress Raw Reads

Compress all .fastq files to .fastq.gz to save disk space.
Wait for the terminal prompt to return before proceeding — this takes several minutes.

```bash
gzip ~/Documents/GitHub/myProject/myProject/raw_reads/*.fastq
```

Verify compression completed successfully:
```bash
ls ~/Documents/GitHub/myProject/myProject/raw_reads/
# Expected: 16 files total (2 per sample x 8 samples), all ending in .fastq.gz
```

---

## Step 4: Quality Control with FastQC

FastQC analyzes raw sequencing reads and generates an HTML report flagging potential
issues such as low base quality scores, adapter contamination, GC content anomalies,
and sequence duplication levels.

```bash
fastqc ~/Documents/GitHub/myProject/myProject/raw_reads/*.fastq.gz \
    --outdir ~/Documents/GitHub/myProject/myProject/qc_results/ \
    --threads 4
```

Open the results folder to view reports in Finder:
```bash
open ~/Documents/GitHub/myProject/myProject/qc_results/
```
Double-click any `.html` file to view the full QC report in your browser.

---

## Sample Summary Table

| # | Accession | Lineage | Resistance Profile | Why This Sample Was Chosen |
|---|---|---|---|---|
| 1 | ERR4796930 | TBD | TBD | One of three initial CRyPTIC dataset samples included to establish a baseline representation of the overall dataset and to test the QC pipeline before scaling up to additional samples |
| 2 | ERR4796931 | TBD | TBD | One of three initial CRyPTIC dataset samples included to establish a baseline representation of the overall dataset and to test the QC pipeline before scaling up to additional samples |
| 3 | ERR4796932 | TBD | TBD | One of three initial CRyPTIC dataset samples included to establish a baseline representation of the overall dataset and to test the QC pipeline before scaling up to additional samples |
| 4 | ERR4796935 | Lineage 4 (Euro-American) | TBD | Lineage 4 is the most globally distributed *M. tuberculosis* lineage, responsible for the majority of TB cases in Europe, the Americas, and parts of Africa. Including this sample ensures the most epidemiologically prevalent lineage is represented in the phylogeny and provides a well-studied backbone for the tree |
| 5 | ERR4796940 | TBD | Drug-susceptible | A fully drug-susceptible isolate is essential as both a phylogenetic and phenotypic baseline. On the tree it represents the ancestral resistance state, allowing us to map the direction and order in which resistance mutations arose across all other samples in the dataset |
| 6 | ERR4796945 | Lineage 2 (Beijing) | MDR-associated | Lineage 2 (Beijing) is one of the most clinically significant lineages due to its disproportionate association with the global emergence and spread of MDR-TB. This isolate was chosen to represent a lineage under intense antibiotic selective pressure and is critical for studying the convergent evolution of resistance across geographically distinct populations |
| 7 | ERR4796950 | Lineage 1 (Indo-Oceanic) | TBD | Lineage 1 is considered the most ancestral of the four major *M. tuberculosis* lineages and is predominantly found in South and Southeast Asia and East Africa. Including this sample is essential for rooting the phylogenetic tree, as it provides the deepest branching point and serves as a reference for interpreting evolutionary relationships among all other isolates |
| 8 | ERR4796960 | TBD | MDR | This isolate satisfies the clinical definition of MDR-TB — resistant to both rifampicin and isoniazid, the two most important first-line TB antibiotics. MDR-TB is a central focus of the CRyPTIC paper, and this sample allows direct genomic comparison of resistance-associated variants against the drug-susceptible baseline (ERR4796940), making it key to any phylogenetic analysis of resistance evolution |

> **Note:** Lineage assignments and resistance profiles marked TBD will be confirmed
> by cross-referencing the CRyPTIC metadata file (CRyPTIC_reuse_table_20211019.csv)
> available via the paper's FTP site in a future step.

---

## QC Results & Observations

*Fill in after reviewing FastQC .html reports in the qc_results/ folder*

| Sample | Per Base Quality | Adapter Content | GC Content | Flags |
|---|---|---|---|---|
| ERR4796930_1 | | | | |
| ERR4796930_2 | | | | |
| ERR4796931_1 | | | | |
| ERR4796931_2 | | | | |
| ERR4796932_1 | | | | |
| ERR4796932_2 | | | | |
| ERR4796935_1 | | | | |
| ERR4796935_2 | | | | |
| ERR4796940_1 | | | | |
| ERR4796940_2 | | | | |
| ERR4796945_1 | | | | |
| ERR4796945_2 | | | | |
| ERR4796950_1 | | | | |
| ERR4796950_2 | | | | |
| ERR4796960_1 | | | | |
| ERR4796960_2 | | | | |

**Expected values for *M. tuberculosis* Illumina data:**
- Per base quality: Phred score ≥ 28 (green in FastQC)
- GC content: ~65% (M. tb has a high GC genome — higher than the human average of ~41%)
- Adapter content: Should be minimal; if flagged, trim with Trimmomatic before proceeding

---

## Next Steps
- [ ] Review all FastQC `.html` reports and fill in QC observations table above
- [ ] Trim adapters if needed using Trimmomatic or Fastp
- [ ] Download CRyPTIC metadata file to confirm lineage and resistance calls for each sample
- [ ] Read mapping to *M. tuberculosis* H37Rv reference genome (NC_000962.3)
- [ ] Variant calling (SNP identification)
- [ ] Multiple sequence alignment
- [ ] Phylogenetic tree construction
- [ ] Push all updates to GitHub
