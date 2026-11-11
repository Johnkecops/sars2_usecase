# Module: 04_models_viz.R
# Purpose: (a) chunked linear model of N-fraction ~ time + platform (biglm idea, slides + kel09 report)
#          (b) overplotting fixes on the PCA scatter (slides 13-17)
# Author: Dr. Arli Aditya Parikesit | Date: 2026
# Parameters: chunk_rows for the stream. Engine = biglm if installed, else normal equations.
source("scripts/helpers.R"); suppressPackageStartupMessages({library(bigmemory); library(ggplot2)})
set.seed(CFG$seed)
bm  <- attach.big.matrix(P(sprintf("features_k%d.desc", CFG$kmer)))
acc <- fread(P("feature_index.tsv"))$accession
meta <- fread(P(CFG$meta_file), select = unname(CFG$cols[c("accession", "date", "platform")]))
setnames(meta, c("accession", "date", "platform"))
mt <- data.table(accession = acc, n_frac = bm[, 2])[meta, on = "accession"]
# ponytail: this join holds 3 small columns in RAM; beyond ~100M rows do it in a database (DuckDB/SQL).
mt[, t_years := as.numeric(as.Date(date) - as.Date("2020-01-01")) / 365.25]
mt_file <- P("model_table.tsv"); fwrite(mt[, .(n_frac, t_years, platform)], mt_file, sep = "\t")

PLAT <- c("Illumina", "IonTorrent", "Nanopore")
prep <- function(d) { d$platform <- factor(d$platform, levels = PLAT); as.data.frame(d) }
fit  <- stream_fit(mt_file, n_frac ~ t_years + platform, prep, chunk_rows = 5000)
print(fit)

# validation on this small set: full in-RAM lm must agree with the streamed fit
ref <- coef(lm(n_frac ~ t_years + platform, prep(copy(mt))))
cat("max |streamed - lm| coefficient difference:", signif(max(abs(fit$estimate - ref)), 3), "\n")
fwrite(fit, "results/model_coefficients.csv")

# overplotting: full scatter vs small dots vs 5% sample vs hexbin
p <- readRDS("results/pca_all.rds"); smp <- p[sample(.N, ceiling(.N * 0.05))]
g <- function(d, geom, title) ggplot(d, aes(PC1, PC2)) + geom + labs(title = title) + theme_minimal(base_size = 11)
ggsave("results/figures/overplot_full.png",   g(p, geom_point(alpha = 1, size = 1.5), sprintf("All points (n=%d)", nrow(p))), width = 5, height = 4, dpi = 120)
ggsave("results/figures/overplot_small.png",  g(p, geom_point(alpha = .15, size = .3), "Small, transparent dots"), width = 5, height = 4, dpi = 120)
ggsave("results/figures/overplot_sample.png", g(smp, geom_point(alpha = .6, size = .8), sprintf("5%% sample (n=%d)", nrow(smp))), width = 5, height = 4, dpi = 120)
ggsave("results/figures/overplot_hex.png",    g(p, geom_hex(bins = 40), "Hexbin density"), width = 5, height = 4, dpi = 120)
cat("figures written to results/figures\n")
