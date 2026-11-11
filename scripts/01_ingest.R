# Module: 01_ingest.R
# Purpose: compare ways to load the metadata table (slides 5, 10-12, 25-26) and record time + R heap change
# Author: Dr. Arli Aditya Parikesit | Date: 2026
# Parameters: CFG (helpers.R). select= keeps only the columns we need (the cheapest memory saving).
# Note: ff is called with ff:: (not library) because it masks Biostrings functions.
source("scripts/helpers.R")
f <- P(CFG$meta_file); keep <- unname(CFG$cols[c("accession", "date", "country", "lineage")])

r <- list(
  measure("base read.delim (all columns)",  read.delim(f, stringsAsFactors = FALSE)),
  measure("data.table fread (all columns)", fread(f)),
  measure("data.table fread (select 4 cols)", fread(f, select = keep)),
  measure("ff read.table.ffdf (on disk)",   ff::read.table.ffdf(file = f, header = TRUE, sep = "\t", quote = ""))
)
res <- rbindlist(lapply(r, `[[`, "row"))
res[, file_mb := round(file.size(f) / 1e6, 2)]
dir.create("results", showWarnings = FALSE); fwrite(res, "results/ingest_benchmark.csv")
print(res)
cat("R", as.character(getRversion()), "| data.table", as.character(packageVersion("data.table")),
    "| ff", as.character(packageVersion("ff")), "\n")
