# Module: 02_sample.R
# Purpose: how wrong is a subsample? Compare lineage proportions: full vs simple vs stratified sample (slides 6-8, 16)
# Author: Dr. Arli Aditya Parikesit | Date: 2026
# Parameters: frac (sampling fraction), seed (CFG$seed)
source("scripts/helpers.R"); suppressPackageStartupMessages(library(dplyr))
set.seed(CFG$seed); frac <- 0.05
d <- fread(P(CFG$meta_file), select = unname(CFG$cols[c("accession", "date", "country", "lineage")]))
setnames(d, c("acc", "date", "country", "lineage"))
d[, `:=`(month = to_month(date), country = to_country(country))]

simple     <- as.data.table(sample_frac(as.data.frame(d), frac))                    # dplyr, as in the deck
stratified <- d[, .SD[sample(.N, max(1L, ceiling(.N * frac)))], by = .(month, country)]
prop <- function(x) x[, .N, by = lineage][, p := N / sum(N)][order(lineage)][, .(lineage, p)]
cmp <- Reduce(function(a, b) merge(a, b, by = "lineage"),
              list(setnames(prop(d), "p", "full"), setnames(prop(simple), "p", "simple"),
                   setnames(prop(stratified), "p", "stratified")))
cmp[, `:=`(err_simple = abs(simple - full), err_strat = abs(stratified - full))]
print(cmp[, lapply(.SD, round, 4), .SDcols = is.numeric, by = lineage])
cat(sprintf("rows: full %d | simple %d | stratified %d | max abs error simple %.4f, stratified %.4f\n",
            nrow(d), nrow(simple), nrow(stratified), max(cmp$err_simple), max(cmp$err_strat)))
cat("Caution: rare lineages (low N) are where a subsample misleads most.\n")
fwrite(cmp, "results/sampling_error.csv")
