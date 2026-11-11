# Scalable SARS-CoV-2 genome analysis in R (BI2106 Session 12 use case)

Author: Dr. Arli Aditya Parikesit, i3L University. Companion to the deck
`Session 12 Use Case - SARS-CoV-2 Large Genomes.pptx`.

A metadata table plus genome FASTA that is too big for `read.csv()` and `lm()` is processed with sampling,
chunked streaming, a file-backed `big.matrix`, a streamed linear model, and a Shiny dashboard.
The bundled data are SYNTHETIC (20,000 genomes of 1,500 bp, simulated lineages and N-fractions).
Conclusions from them apply to the code and memory behaviour, not to SARS-CoV-2 biology.

## Run order in RStudio
1. Open RStudio, set the working directory to this folder (Session > Set Working Directory).
2. Install packages once:
   `install.packages(c("data.table","bigmemory","ff","shiny","bslib","ggplot2","dplyr","hexbin","DT","lobstr","BiocManager"))`
   then `BiocManager::install("Biostrings")`. Optional, used if present: `pryr`, `biglm`.
3. Run `scripts/00_make_synthetic.R` (skip if you use real data), then `01` to `04` in order.
4. Open `app.R` and press Run App.

`data/genomes.fasta` and `data/features_k3.bin` are not in git (large, regenerable): run `scripts/00` and `scripts/03` to create them.

Every script starts with `source("scripts/helpers.R")`, so run them from this folder.
Do not name the script folder `R/`: Shiny auto-sources everything in an `R/` folder.

## Using real NCBI data
Download SARS-CoV-2 genomes and metadata from NCBI Virus or the `datasets` tool; see the NCBI how-to
https://www.ncbi.nlm.nih.gov/datasets/docs/how-tos/virus/get-sars2-genomes/ . Save as `data/metadata.tsv`
(tab-separated) and `data/genomes.fasta`. Then edit `CFG$cols` in `scripts/helpers.R`.
The column names there are my best guess, not verified against a current download. Check your header first.
The metadata `Accession` must match the FASTA record IDs, and the FASTA order defines the matrix row order.
Set `SARS2_DATA_DIR` to use another data folder.

## What was tested, and what was not
Tested (Linux sandbox, R 4.3.3, data.table 1.14.10, bigmemory 4.6.1, ff 4.0.12, Biostrings 2.70.2, shiny 1.8.0,
ggplot2 3.4.4, dplyr 1.1.4): scripts 00 to 04 and every dashboard tab on the synthetic data.
Streamed model matches `lm()` to 5.3e-16 on this data.

NOT tested: the `biglm` code path (CRAN was unreachable in the sandbox, so the normal-equation fallback ran),
`pryr` (lobstr supplied the heap numbers), `ffbase`, real NCBI data, and timings at full scale on a laptop.
Timings on the 20,000-row demo are milliseconds and not informative; compare heap change instead.

## Known limits
- k-mer size is 3 in the dashboard; set `CFG$kmer` and rerun script 03 to change it (the app reads k = 3 files).
- The normal-equation fallback loses precision with near-collinear predictors; `biglm` is preferred.
- The model-table join keeps three small columns in RAM; very large tables need a database step.
- PC1 and PC2 explain about 13% and 7% of k-mer variance in the synthetic set. No lineage separation is claimed.

AI Assistance Disclaimer: This codebase was developed with the assistance of Claude Code. While the AI provided code generation, debugging, and structural support, the human developer maintains full responsibility for reviewing, testing, and maintaining all content and functionality
