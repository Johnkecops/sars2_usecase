# Module: 00_make_synthetic.R
# Purpose: build a SYNTHETIC SARS-CoV-2-like collection (metadata.tsv + genomes.fasta)
#          so the whole use case runs offline. NOT real biology.
# Author: Dr. Arli Aditya Parikesit | Date: 2026
# Parameters: CFG$n_genomes, CFG$genome_len, CFG$seed (see helpers.R)
# Built-in "ground truth" (for checking the model in 04_models_viz.R):
#   N-fraction mean = platform base rate * exp(-0.012 * months since 2020-01)
#   base rates: Illumina 0.004, IonTorrent 0.008, Nanopore 0.010
#   -> the fitted time slope should be negative and Nanopore > Illumina.
# Real data: skip this script, download from NCBI Virus (see README) and edit CFG$cols.
source("scripts/helpers.R")
suppressPackageStartupMessages(library(Biostrings))
set.seed(CFG$seed); N <- CFG$n_genomes; L <- CFG$genome_len
dir.create(CFG$data_dir, showWarnings = FALSE)

lineages  <- paste0("Lineage_", LETTERS[1:6]); emerge <- c(0, 6, 14, 22, 30, 40)
countries <- c("Indonesia", "USA", "India", "Brazil", "United Kingdom", "South Africa", "Japan")
month <- sample(0:53, N, TRUE, prob = dnorm(0:53, 20, 14))
W <- sapply(emerge, function(e) dnorm(month, e + 6, 6)) + 1e-9
lin <- pmin(1L + rowSums(runif(N) > t(apply(W / rowSums(W), 1, cumsum))), 6L)
plat <- ifelse(runif(N) < 0.3 + 0.004 * month, "Nanopore",
               sample(c("Illumina", "IonTorrent"), N, TRUE, c(.75, .25)))
mu <- c(Illumina = 0.004, Nanopore = 0.010, IonTorrent = 0.008)[plat] * exp(-0.012 * month)
nfrac <- pmin(rexp(N, 1 / mu), 0.2)

# SARS-CoV-2-like base composition (approximate), nudged per lineage so k-mers carry signal
base <- c(A = .30, C = .18, G = .20, T = .32)
pl <- t(sapply(1:6, function(k) { p <- pmax(base + rnorm(4, 0, .015), .05); p / sum(p) }))
acc <- sprintf("SYN%06d", seq_len(N))

meta <- data.table(acc, as.character(as.Date("2020-01-01") + round(month * 30.44) + sample(0:29, N, TRUE)),
                   paste0(sample(countries, N, TRUE), ": region"), lineages[lin], L, plat)
setnames(meta, unname(CFG$cols[c("accession", "date", "country", "lineage", "length", "platform")]))
fwrite(meta, P(CFG$meta_file), sep = "\t")

fa <- P(CFG$fasta_file); if (file.exists(fa)) file.remove(fa)
for (s in seq(1L, N, by = CFG$chunk_genomes)) {
  i <- s:min(s + CFG$chunk_genomes - 1L, N)
  sq <- vapply(i, function(j) {
    v <- c("A", "C", "G", "T")[sample.int(4L, L, TRUE, pl[lin[j], ])]
    k <- round(nfrac[j] * L); if (k > 0) v[sample.int(L, k)] <- "N"
    paste(v, collapse = "")
  }, "")
  x <- DNAStringSet(sq); names(x) <- acc[i]
  writeXStringSet(x, fa, append = s > 1L)
}
cat("wrote", N, "genomes of", L, "bp:", P(CFG$meta_file), "and", fa, "\n")
