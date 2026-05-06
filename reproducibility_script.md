# M. tuberculosis Phylogenetic Analysis: Reproducibility Script
**Aidan Guyse | May 2026**

---

## Software Cheatsheet

| Software | Description | Strengths | Weaknesses | Assumptions | User Choices |
|---|---|---|---|---|---|
| **FastQC v0.12** | QC for raw Illumina reads | Fast, visual output | Does not trim reads | Illumina short reads | Phred >30, GC ~65% |
| **BWA MEM v0.7.17** | Maps paired-end reads to reference | Fast, standard for M. tb WGS | Requires reference genome | Illumina paired-end reads | 4 threads, H37Rv reference |
| **SAMtools v1.17** | Processes BAM/SAM files | Efficient, widely supported | Command line only | Sorted, indexed BAM input | Most frequent base called |
| **MAFFT v7.526** | Multiple sequence alignment | Handles very long sequences | FFT-NS-2 less accurate than iterative methods | Nucleotide input | `--retree 2`, 4 threads |
| **IQ-TREE v2** | Maximum likelihood inference | Fast UFBoot, wide model support | UFBoot can overestimate support | Aligned sequences, user specified model | GTR+G, 1000 UFBoot, ERR4796950 outgroup |
| **RAxML-NG v1.2** | Maximum likelihood inference | Standard bootstrap, independent cross check | Slower than IQ-TREE | Aligned sequences, user specified model | GTR+G, 1000 bootstrap, ERR4796950 outgroup |
| **MrBayes v3.2.7** | Bayesian phylogenetic inference | Posterior probabilities | Very slow for large datasets | MCMC convergence required | 1,000,000 generations, 4 chains, 25% burnin |
| **ASTRAL-III v5.15** | Coalescent species tree inference | Handles ILS, fast | Invalid for clonal organisms | Sexual reproduction, recombination | Demonstrated on mammal toy dataset only |
| **ape v5.7** | Phylogenetic analysis in R | NJ, tree manipulation, visualization | Limited ML capabilities | Rooted tree for plotting | TN93 for NJ, ERR4796950 outgroup |
| **phangorn v2.11** | Phylogenetic analysis in R | Parsimony inference | Slower than dedicated ML tools | phyDat format input | NNI branch swapping, NJ starting tree |

---

## To Do Before Starting

- [ ] Install all software listed in the cheatsheet above
- [ ] Download raw reads for all 29 isolates from NCBI SRA, accessions are listed in Table 1 of the paper
- [ ] Download H37Rv reference genome accession GCF_000195955.2 and save it to `reference/`
- [ ] Clone the repo: `git clone https://github.com/aguyse/myProject.git`
- [ ] Install R packages: `install.packages(c("ape", "phangorn"))`

---

## Step 1: Set Up Directory Structure

```bash
mkdir -p ~/Documents/GitHub/myProject/myProject/data/bam
mkdir -p ~/Documents/GitHub/myProject/myProject/data/consensus
mkdir -p ~/Documents/GitHub/myProject/myProject/data/alignment
mkdir -p ~/Documents/GitHub/myProject/myProject/data/ml_tree
mkdir -p ~/Documents/GitHub/myProject/myProject/data/raxml_tree
mkdir -p ~/Documents/GitHub/myProject/myProject/data/bayesian
mkdir -p ~/Documents/GitHub/myProject/myProject/data/coalescent
mkdir -p ~/Documents/GitHub/myProject/myProject/qc_results
mkdir -p ~/Documents/GitHub/myProject/myProject/scripts

cd ~/Documents/GitHub/myProject
git init
git remote add origin https://github.com/aguyse/myProject.git
```

---

## Step 2: Quality Control

FastQC on all raw reads. All 29 isolates should pass Phred >30 and GC content around 65%.

```bash
fastqc raw_reads/*.fastq.gz -o ~/Documents/GitHub/myProject/myProject/qc_results/
```

```bash
cd ~/Documents/GitHub/myProject
git add .
git commit -m "Add FastQC quality control reports"
git push
```

---

## Step 3: Index Reference Genome

Only needs to be done once before mapping.

```bash
bwa index ~/Documents/GitHub/myProject/myProject/reference/H37Rv.fasta
```

---

## Step 4: Map Reads to Reference Using BWA MEM

BWA MEM was chosen over Bowtie2 as it is the standard aligner for M. tb WGS pipelines. Each sample is mapped, sorted, and indexed.

```bash
for sample in ERR4812134 ERR4828611 ERR4812438 ERR8699122 ERR4831359 \
              ERR8975594 ERR4809142 ERR2510615 ERR4813066 ERR4828858 \
              ERR4829143 ERR4813873 ERR4796544 ERR4830996 ERR4798070 \
              ERR4796794 ERR4831425 ERR8975515 ERR2510399 ERR8665736 \
              ERR8665737 ERR4796930 ERR4796931 ERR4796932 ERR4796935 \
              ERR4796940 ERR4796945 ERR4796950 ERR4796960; do

    bwa mem -t 4 \
        ~/Documents/GitHub/myProject/myProject/reference/H37Rv.fasta \
        raw_reads/${sample}_1.fastq.gz \
        raw_reads/${sample}_2.fastq.gz | \
    samtools sort -o ~/Documents/GitHub/myProject/myProject/data/bam/${sample}.bam

    samtools index ~/Documents/GitHub/myProject/myProject/data/bam/${sample}.bam
    echo "Done: $sample"
done

# Verify all 29 BAMs were created
ls ~/Documents/GitHub/myProject/myProject/data/bam/*.bam | wc -l
# Expected output: 29
```

```bash
cd ~/Documents/GitHub/myProject
git add .
git commit -m "Add BWA mapped BAM files for all 29 isolates"
git push
```

---

## Step 5: Generate Consensus Sequences

SAMtools consensus calls the most frequent base at each position for each sample.

```bash
for sample in ERR4812134 ERR4828611 ERR4812438 ERR8699122 ERR4831359 \
              ERR8975594 ERR4809142 ERR2510615 ERR4813066 ERR4828858 \
              ERR4829143 ERR4813873 ERR4796544 ERR4830996 ERR4798070 \
              ERR4796794 ERR4831425 ERR8975515 ERR2510399 ERR8665736 \
              ERR8665737 ERR4796930 ERR4796931 ERR4796932 ERR4796935 \
              ERR4796940 ERR4796945 ERR4796950 ERR4796960; do

    samtools consensus -f fasta \
        ~/Documents/GitHub/myProject/myProject/data/bam/${sample}.bam \
        -o ~/Documents/GitHub/myProject/myProject/data/consensus/${sample}.fasta

    sed -i '' "s/>.*/>$sample/" \
        ~/Documents/GitHub/myProject/myProject/data/consensus/${sample}.fasta
done

# Combine all 29 into one multi-FASTA
cat ~/Documents/GitHub/myProject/myProject/data/consensus/*.fasta > \
    ~/Documents/GitHub/myProject/myProject/data/alignment/all_samples_combined.fasta

# Verify all 29 are present
grep ">" ~/Documents/GitHub/myProject/myProject/data/alignment/all_samples_combined.fasta | wc -l
# Expected output: 29
```

```bash
cd ~/Documents/GitHub/myProject
git add .
git commit -m "Add consensus FASTA sequences and combined multi-FASTA"
git push
```

---

## Step 6: Multiple Sequence Alignment Using MAFFT

MUSCLE was attempted first but crashed because whole genomes at 4.4 Mb exceed its length limit. MAFFT FFT-NS-2 handles arbitrarily long sequences without a length ceiling. Any negative distance warnings during the run are normal for highly similar sequences and do not affect the output.

```bash
mafft --retree 2 --nuc --thread 4 \
    ~/Documents/GitHub/myProject/myProject/data/alignment/all_samples_combined.fasta \
    > ~/Documents/GitHub/myProject/myProject/data/alignment/mafft_aligned.fasta

# Verify 29 sequences in the alignment
grep ">" ~/Documents/GitHub/myProject/myProject/data/alignment/mafft_aligned.fasta | wc -l
# Expected output: 29
```

```bash
cd ~/Documents/GitHub/myProject
git add .
git commit -m "Add MAFFT whole genome alignment"
git push
```

---

## Step 7: Neighbor Joining and Maximum Parsimony Trees

Both trees are built in R. NJ uses TN93 distances which is appropriate for M. tb given its 65% GC content since simpler models like JC69 assume equal base frequencies, which is violated here. MP uses the NJ tree as a starting topology with NNI branch swapping.

```r
library(ape)
library(phangorn)

dna <- read.dna("~/Documents/GitHub/myProject/myProject/data/alignment/mafft_aligned.fasta",
                format = "fasta")
dna_phyDat <- as.phyDat(dna)

# NJ tree
dist_matrix <- dist.dna(dna, model = "TN93", pairwise.deletion = TRUE)
nj_tree <- nj(dist_matrix)
nj_tree_rooted <- root(nj_tree, outgroup = "ERR4796950", resolve.root = TRUE)
write.tree(nj_tree_rooted,
    "~/Documents/GitHub/myProject/myProject/data/nj_tree.nwk")

# MP tree
# Expected output: "Final p-score 11938 after 6 nni operations"
mp_tree <- optim.parsimony(nj_tree, dna_phyDat)
mp_tree_rooted <- root(mp_tree, outgroup = "ERR4796950", resolve.root = TRUE)
write.tree(mp_tree_rooted,
    "~/Documents/GitHub/myProject/myProject/data/mp_tree.nwk")
```

```bash
cd ~/Documents/GitHub/myProject
git add .
git commit -m "Add NJ and MP phylogenetic trees"
git push
```

---

## Step 8: Maximum Likelihood Trees

Two independent ML tools are run as a cross check. Both use GTR+G. GTR is chosen because M. tb violates equal base frequency assumptions of simpler models, and the +G parameter accounts for rate variation across sites. If both tools arrive at the same tree the result is much more trustworthy.

**IQ-TREE, takes about 38 seconds:**
```bash
# Expected output: lnL = -6,094,539.401
iqtree -s ~/Documents/GitHub/myProject/myProject/data/alignment/mafft_aligned.fasta \
    -o ERR4796950 \
    -m GTR+G \
    -B 1000 \
    -T 4 \
    --prefix ~/Documents/GitHub/myProject/myProject/data/ml_tree/mtb_ml \
    --redo
```

**RAxML-NG, takes about 2 hours:**
```bash
# Expected output: lnL = -6,094,537.794
raxml-ng --msa ~/Documents/GitHub/myProject/myProject/data/alignment/mafft_aligned.fasta \
    --model GTR+G \
    --outgroup ERR4796950 \
    --bs-trees 1000 \
    --all \
    --threads 4 \
    --redo \
    --prefix ~/Documents/GitHub/myProject/myProject/data/raxml_tree/mtb_raxml
```

```bash
cd ~/Documents/GitHub/myProject
git add .
git commit -m "Add IQ-TREE and RAxML-NG maximum likelihood trees"
git push
```

---

## Step 9: Bayesian Inference Using MrBayes

Watch the ASDSF column as it runs. It needs to drop below 0.01 to indicate convergence. In this analysis it reached 0.000000 well before the run completed. Expect this to take roughly 2 hours.

First convert the alignment to NEXUS format in R:

```r
library(ape)
dna <- read.dna("~/Documents/GitHub/myProject/myProject/data/alignment/mafft_aligned.fasta",
                format = "fasta")
write.nexus.data(dna,
    file = "/Users/aidanguyse/Documents/GitHub/myProject/myProject/data/bayesian/mtb_alignment.nex",
    format = "dna")
```

Then create the MrBayes input file and run:

```bash
cat > /Users/aidanguyse/Documents/GitHub/myProject/myProject/data/bayesian/mrbayes_run.nex << 'EOF'
begin mrbayes;
  set autoclose=yes;
  execute /Users/aidanguyse/Documents/GitHub/myProject/myProject/data/bayesian/mtb_alignment.nex;
  outgroup ERR4796950;
  lset nst=6 rates=gamma;
  mcmc ngen=1000000 samplefreq=1000 printfreq=10000 nchains=4;
  sump burnin=250;
  sumt burnin=250;
end;
EOF

mb /Users/aidanguyse/Documents/GitHub/myProject/myProject/data/bayesian/mrbayes_run.nex
```

```bash
cd ~/Documents/GitHub/myProject
git add .
git commit -m "Add MrBayes Bayesian inference output and consensus tree"
git push
```

---

## Step 10: ASTRAL on Toy Dataset

ASTRAL is not applicable to M. tb because it is strictly clonal and all genomic regions share one evolutionary history, so gene tree discordance is absent. It is demonstrated here on a published mammalian dataset to illustrate how the method works.

```bash
astral -i ~/miniconda3/pkgs/astral-tree-*/share/astral*/test_data/song_mammals.424.gene.tre \
    -o ~/Documents/GitHub/myProject/myProject/data/coalescent/astral_mammals.tre \
    2>&1 | tee ~/Documents/GitHub/myProject/myProject/data/coalescent/astral_mammals.log
```

```bash
cd ~/Documents/GitHub/myProject
git add .
git commit -m "Add ASTRAL coalescent demonstration on mammal toy dataset"
git push
```

---

## Step 11: Tree Visualization

All trees are annotated by resistance category. Note that ERR8665736 and ERR8665737 are missing from the metadata CSV and must be hardcoded manually as INH Monoresistant and Pre-XDR respectively. Without this fix they will appear black on the tree.

```r
library(ape)

metadata <- read.csv("~/Documents/GitHub/myProject/myProject/data/isolate_metadata.csv")

color_map <- c(
  "Drug Susceptible"  = "#4DAF4A",
  "INH Monoresistant" = "#BFAF3A",
  "RIF Monoresistant" = "#E07B39",
  "MDR"               = "#9B59B6",
  "Pre-XDR"           = "#C0392B",
  "XDR"               = "#7B0000"
)

# Helper to get tip colors with manual fixes for the two missing isolates
get_colors <- function(tree) {
  tip_colors <- metadata$RESISTANCE_CATEGORY[match(tree$tip.label, metadata$ENA_RUN)]
  tip_colors[tree$tip.label == "ERR8665736"] <- "INH Monoresistant"
  tip_colors[tree$tip.label == "ERR8665737"] <- "Pre-XDR"
  color_map[tip_colors]
}

# NJ tree
nj_tree <- read.tree("~/Documents/GitHub/myProject/myProject/data/nj_tree.nwk")
pdf("~/Documents/GitHub/myProject/myProject/qc_results/nj_tree_annotated.pdf", width = 18, height = 14)
par(mar = c(2, 2, 4, 12), xpd = TRUE)
plot(nj_tree, tip.color = get_colors(nj_tree),
     main = "Neighbor Joining Tree - M. tuberculosis (TN93)", cex = 0.9,
     x.lim = c(0, max(node.depth.edgelength(nj_tree)) * 1.5))
legend(x = par("usr")[2] + 0.000001, y = par("usr")[4],
       legend = names(color_map), col = color_map,
       pch = 19, cex = 0.85, title = "Resistance Category", bty = "n")
dev.off()

# MP tree uses cladogram type because parsimony produces no branch lengths
mp_tree <- read.tree("~/Documents/GitHub/myProject/myProject/data/mp_tree.nwk")
pdf("~/Documents/GitHub/myProject/myProject/qc_results/mp_tree_annotated.pdf", width = 18, height = 14)
par(mar = c(2, 2, 4, 12), xpd = TRUE)
plot(mp_tree, type = "cladogram", tip.color = get_colors(mp_tree),
     main = "Maximum Parsimony Tree - M. tuberculosis", cex = 0.9)
legend(x = par("usr")[2] + 0.01, y = par("usr")[4],
       legend = names(color_map), col = color_map,
       pch = 19, cex = 0.85, title = "Resistance Category", bty = "n")
dev.off()

# IQ-TREE ML tree with bootstrap values shown at nodes
ml_tree <- read.tree("~/Documents/GitHub/myProject/myProject/data/ml_tree/mtb_ml.treefile")
ml_tree_rooted <- root(ml_tree, outgroup = "ERR4796950", resolve.root = TRUE)
pdf("~/Documents/GitHub/myProject/myProject/qc_results/ml_tree_annotated.pdf", width = 18, height = 14)
par(mar = c(2, 2, 4, 12), xpd = TRUE)
plot(ml_tree_rooted, tip.color = get_colors(ml_tree_rooted),
     main = "IQ-TREE Maximum Likelihood Tree - M. tuberculosis (GTR+G)", cex = 0.9,
     x.lim = c(0, max(node.depth.edgelength(ml_tree_rooted)) * 1.5))
nodelabels(ml_tree_rooted$node.label, adj = c(1.1, -0.3), frame = "none", cex = 0.6)
legend(x = par("usr")[2] + 0.000001, y = par("usr")[4],
       legend = names(color_map), col = color_map,
       pch = 19, cex = 0.85, title = "Resistance Category", bty = "n")
dev.off()

# RAxML tree with bootstrap values shown at nodes
raxml_tree <- read.tree("~/Documents/GitHub/myProject/myProject/data/raxml_tree/mtb_raxml.raxml.support")
raxml_tree_rooted <- root(raxml_tree, outgroup = "ERR4796950", resolve.root = TRUE)
pdf("~/Documents/GitHub/myProject/myProject/qc_results/raxml_tree_annotated.pdf", width = 18, height = 14)
par(mar = c(2, 2, 4, 12), xpd = TRUE)
plot(raxml_tree_rooted, tip.color = get_colors(raxml_tree_rooted),
     main = "RAxML-NG Maximum Likelihood Tree - M. tuberculosis (GTR+G)", cex = 0.9,
     x.lim = c(0, max(node.depth.edgelength(raxml_tree_rooted)) * 1.5))
nodelabels(raxml_tree_rooted$node.label, adj = c(1.1, -0.3), frame = "none", cex = 0.6)
legend(x = par("usr")[2] + 0.000001, y = par("usr")[4],
       legend = names(color_map), col = color_map,
       pch = 19, cex = 0.85, title = "Resistance Category", bty = "n")
dev.off()

# MrBayes tree
bayes_tree <- read.nexus("~/Documents/GitHub/myProject/myProject/data/bayesian/mtb_alignment.nex.con.tre")
bayes_tree_rooted <- root(bayes_tree, outgroup = "ERR4796950", resolve.root = TRUE)
pdf("~/Documents/GitHub/myProject/myProject/qc_results/bayes_tree_annotated.pdf", width = 18, height = 14)
par(mar = c(2, 2, 4, 12), xpd = TRUE)
plot(bayes_tree_rooted, tip.color = get_colors(bayes_tree_rooted),
     main = "Bayesian Consensus Tree - M. tuberculosis MrBayes GTR+G", cex = 0.9,
     x.lim = c(0, max(node.depth.edgelength(bayes_tree_rooted)) * 1.5))
legend(x = par("usr")[2] + 0.000001, y = par("usr")[4],
       legend = names(color_map), col = color_map,
       pch = 19, cex = 0.85, title = "Resistance Category", bty = "n")
dev.off()

# Cophylogeny comparing IQ-TREE vs RAxML
# Line crossings reflect arbitrary tip display order, not topological disagreement
# Both trees have near identical log-likelihoods confirming the same optimal topology
pdf("~/Documents/GitHub/myProject/myProject/qc_results/cophylogeny_iqtree_raxml.pdf",
    width = 20, height = 14)
association <- cbind(ml_tree_rooted$tip.label, ml_tree_rooted$tip.label)
cophyloplot(ml_tree_rooted, raxml_tree_rooted, assoc = association,
            length.line = 4, space = 28, gap = 3, col = "steelblue", lwd = 0.8,
            main = "Cophylogenetic Comparison: IQ-TREE vs RAxML-NG")
dev.off()

# MrBayes convergence plot, ASDSF should drop to 0 rapidly
log_lines <- readLines("~/Documents/GitHub/myProject/myProject/data/bayesian/mtb_alignment.nex.mcmc")
asdsf_lines <- grep("Average standard deviation", log_lines, value = TRUE)
asdsf_values <- as.numeric(gsub(".*: ", "", asdsf_lines))
generations <- seq(10000, length.out = length(asdsf_values), by = 10000)
pdf("~/Documents/GitHub/myProject/myProject/qc_results/mrbayes_convergence.pdf", width = 10, height = 6)
plot(generations, asdsf_values, type = "l", lwd = 2, col = "steelblue",
     xlab = "Generation", ylab = "ASDSF",
     main = "MrBayes MCMC Convergence")
abline(h = 0.01, col = "red", lty = 2, lwd = 1.5)
legend("topright", legend = c("ASDSF", "Threshold 0.01"),
       col = c("steelblue", "red"), lty = c(1,2), lwd = 2, bty = "n")
dev.off()
```

---

## Step 12: Commit Everything to GitHub

```bash
cd ~/Documents/GitHub/myProject
git add .
git commit -m "Add all annotated tree visualizations and convergence plot"
git push
```
