# Genomic Analysis of Mycobacterium tuberculosis (CRyPTIC Dataset)

## Project Overview
This project analyzes raw next-generation sequencing (NGS) reads from clinical isolates of the Mycobacterium tuberculosis species complex (MTBC). Data are obtained from the CRyPTIC compendium, which links whole-genome sequencing data to phenotypic resistance profiles for 13 anti-tuberculosis drugs.

### Study Reference
CRyPTIC Consortium (2022)
PMID: 35944069

---

## Dataset Description
- Organism: Mycobacterium tuberculosis complex
- Data type: Paired-end whole genome sequencing
- Source: European Nucleotide Archive (ENA)
- Phenotype data: Resistance profiles for 13 antibiotics

---

## Current Dataset Stage
Raw FASTQ reads

Next Step: Perform quality control using FastQC

---

## Suggested Next Actions
- Run FastQC on all raw FASTQ files and save reports to `qc_results/`
- Aggregate FastQC reports with MultiQC for overview
- Trim adapters and low-quality bases with `fastp` or `TrimGalore` and write trimmed reads to `trimmed_reads/`
- Re-run FastQC on trimmed reads to confirm improvements
- Proceed to read mapping, variant calling, and phenotype association analyses
