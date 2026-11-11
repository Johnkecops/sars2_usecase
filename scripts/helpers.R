# Module: helpers.R
# Purpose: shared config + small helpers for the SARS-CoV-2 large-dataset use case
# Author: Dr. Arli Aditya Parikesit (course BI2106, i3L University)
# Date: 2026
# Parameters:
#   CFG$n_genomes, CFG$genome_len: size of the SYNTHETIC demo set (change to scale up)
#   CFG$cols: map of logical names -> column names in YOUR metadata file.
#             Verify these against the header of your NCBI Virus download.
# Run from the project root (setwd to the folder that holds app.R).

suppressPackageStartupMessages(library(data.table))

CFG <- list(
  data_dir      = Sys.getenv("SARS2_DATA_DIR", "data"),
  meta_file     = "metadata.tsv",
  fasta_file    = "genomes.fasta",
  n_genomes     = 20000L,
  genome_len    = 1500L,     # real SARS-CoV-2 is ~30 kb; kept short so the demo is small
  seed          = 2026L,
  chunk_genomes = 2000L,
  kmer          = 3L,
  cols = c(accession = "Accession", date = "Collection_Date", country = "Geo_Location",
           lineage = "Pangolin", length = "Length", platform = "Sequencing_Technology")
)
P <- function(...) file.path(CFG$data_dir, ...)

# R heap in bytes. pryr::mem_used() is what the slides use; lobstr is the fallback.
mem_used <- function() {
  if (requireNamespace("pryr", quietly = TRUE)) as.numeric(pryr::mem_used())
  else as.numeric(lobstr::mem_used())
}

# Time + heap change of one expression. Returns list(row = data.table, value = result).
measure <- function(label, expr) {
  gc(FALSE); m0 <- mem_used(); t0 <- proc.time()[["elapsed"]]
  val <- force(expr)
  list(row = data.table(label = label,
                        seconds = round(proc.time()[["elapsed"]] - t0, 3),
                        mem_change_mb = round((mem_used() - m0) / 1e6, 2)),
       value = val)
}

# "2021-03-15" or "2021-03" -> first day of month; year-only or junk -> NA
to_month <- function(x) {
  x <- as.character(x); out <- as.Date(rep(NA_character_, length(x)))
  ok <- !is.na(x) & nchar(x) >= 7
  out[ok] <- as.Date(paste0(substr(x[ok], 1, 7), "-01"), format = "%Y-%m-%d")
  out
}
# "USA: California" -> "USA"
to_country <- function(x) sub(":.*$", "", trimws(as.character(x)))

# Chunked least squares: the idea behind biglm (one pass, memory ~ p x p).
# Uses biglm::biglm + update() when biglm is installed, otherwise accumulates
# X'X and X'y (normal equations).
# ponytail: normal equations lose precision when predictors are nearly collinear;
# biglm updates a QR decomposition and is the safer choice. Upgrade path: install biglm.
stream_fit <- function(file, formula, prep, chunk_rows = 1e5, engine = c("auto", "biglm", "normal")) {
  engine <- match.arg(engine)
  if (engine == "auto") engine <- if (requireNamespace("biglm", quietly = TRUE)) "biglm" else "normal"
  hdr <- names(fread(file, nrows = 0)); done <- 0L; fit <- NULL
  XtX <- Xty <- yy <- NULL; n <- 0
  repeat {
    # fread errors when skip is past the last line, so treat that as "no more chunks"
    d <- tryCatch(fread(file, skip = 1L + done, nrows = chunk_rows, header = FALSE, col.names = hdr),
                  error = function(e) NULL)
    if (is.null(d) || !nrow(d)) break
    done <- done + nrow(d); last <- nrow(d) < chunk_rows; d <- prep(d)
    if (engine == "biglm") {
      fit <- if (is.null(fit)) biglm::biglm(formula, d) else update(fit, d)
    } else {
      mf <- model.frame(formula, d); X <- model.matrix(formula, mf); y <- model.response(mf)
      if (is.null(XtX)) { XtX <- crossprod(X); Xty <- crossprod(X, y); yy <- sum(y^2) }
      else { XtX <- XtX + crossprod(X); Xty <- Xty + crossprod(X, y); yy <- yy + sum(y^2) }
      n <- n + nrow(X)
    }
    if (last) break
  }
  if (engine == "biglm") {
    return(data.table(term = names(coef(fit)), estimate = unname(coef(fit)),
                      se = sqrt(diag(vcov(fit))), n = done, engine = engine))
  }
  b <- solve(XtX, Xty); s2 <- (yy - sum(b * Xty)) / (n - length(b))
  data.table(term = rownames(b), estimate = drop(b),
             se = sqrt(diag(solve(XtX)) * s2), n = n, engine = engine)
}
