# Module: 03_features_bigmemory.R
# Purpose: stream the genome FASTA in chunks, compute per-genome features, store them in a
#          file-backed big.matrix (slides 19-24), then PCA fitted on a sample and projected chunk-wise.
# Author: Dr. Arli Aditya Parikesit | Date: 2026
# Parameters: CFG$chunk_genomes, CFG$kmer. Columns: 1 length, 2 n_frac, 3 gc, 4.. k-mer frequencies.
# Biostrings::fasta.index gives random access, so each chunk reads only its own records.
source("scripts/helpers.R")
suppressPackageStartupMessages({library(Biostrings); library(bigmemory)})
fa <- P(CFG$fasta_file); k <- CFG$kmer
idx <- fasta.index(fa); n <- nrow(idx); nk <- 4^k
unlink(file.path(CFG$data_dir, sprintf("features_k%d.%s", k, c("bin", "desc"))))  # allow re-runs
bm <- filebacked.big.matrix(n, 3 + nk, type = "double", backingpath = CFG$data_dir,
                            backingfile = sprintf("features_k%d.bin", k),
                            descriptorfile = sprintf("features_k%d.desc", k))
log <- list(); i <- 0L
for (s in seq(1L, n, by = CFG$chunk_genomes)) {
  rows <- s:min(s + CFG$chunk_genomes - 1L, n)
  m <- measure(sprintf("chunk %d", i <- i + 1L), {
    x  <- readDNAStringSet(idx[rows, ])
    ab <- alphabetFrequency(x)[, c("A", "C", "G", "T", "N")]
    acgt <- rowSums(ab[, c("A", "C", "G", "T")])
    bm[rows, ] <- cbind(width(x), ab[, "N"] / width(x), (ab[, "C"] + ab[, "G"]) / pmax(acgt, 1),
                        oligonucleotideFrequency(x, k, as.prob = TRUE))
  })
  log[[i]] <- m$row
}
fwrite(data.table(accession = sub(" .*", "", idx$desc)), P("feature_index.tsv"))
fwrite(rbindlist(log), "results/feature_chunk_log.csv"); print(rbindlist(log))
cat("big.matrix:", n, "x", ncol(bm), "| backing file MB:", round(file.size(P(sprintf("features_k%d.bin", k))) / 1e6, 1),
    "| in-RAM size if loaded as matrix MB:", round(n * ncol(bm) * 8 / 1e6, 1), "\n")

# PCA: fit on a sample (RAM-safe), project every genome chunk-wise
set.seed(CFG$seed); kc <- 4:(3 + nk)
samp <- sort(sample.int(n, min(5000L, n)))
pc <- prcomp(bm[samp, kc], center = TRUE)
ctr <- pc$center; rot <- pc$rotation[, 1:2]
coords <- do.call(rbind, lapply(split(seq_len(n), ceiling(seq_len(n) / CFG$chunk_genomes)),
                                function(r) sweep(bm[r, kc], 2, ctr) %*% rot))
meta <- fread(P(CFG$meta_file), select = unname(CFG$cols[c("accession", "lineage")])); setnames(meta, c("accession", "lineage"))
out <- data.table(accession = fread(P("feature_index.tsv"))$accession, PC1 = coords[, 1], PC2 = coords[, 2])
out <- merge(out, meta, by = "accession", sort = FALSE)
saveRDS(out, "results/pca_all.rds")
cat("PCA variance explained (sample fit), PC1/PC2:", round(summary(pc)$importance[2, 1:2], 3), "\n")
