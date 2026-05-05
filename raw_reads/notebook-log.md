# Mycobacterium tuberculosis Genomic Analysis
**Name:** Aidan
**Dataset:** CRyPTIC compendium (PubMed: 35944069)
**Description:** Genomic analysis of M. tuberculosis species complex using raw NGS reads to identify drug-resistance mutations for 13 antibiotics.

## Research Log
### Tue Feb 10 23:58:40 CST 2026
- Initialized project repository.
- Manually installed FastQC v0.12.1 due to Homebrew versioning conflicts.
- Performing Quality Control on raw NGS reads.
Command executed: ./FastQC/fastqc raw_reads/*.fastq.gz -o qc_results/
