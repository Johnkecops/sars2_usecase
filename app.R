# Module: app.R
# Purpose: Shiny dashboard to drive the SARS-CoV-2 large-dataset use case from RStudio
#          (Run App button, or shiny::runApp()). Controls: sample fraction, filters,
#          load backend, chunk size, plot mode. A memory/time monitor logs every action.
# Author: Dr. Arli Aditya Parikesit | Date: 2026
# Parameters: set SARS2_DATA_DIR to point to real data; run scripts/01..04 scripts first.
# Needs: shiny, bslib, data.table, dplyr, ggplot2, hexbin, DT, bigmemory, ff (pryr or lobstr)
source("scripts/helpers.R")
suppressPackageStartupMessages({library(shiny); library(bslib); library(ggplot2); library(dplyr)
                                library(DT); library(bigmemory)})

need <- c(P("features_k3.desc"), "results/pca_all.rds", P("model_table.tsv"))
if (!all(file.exists(need))) stop("Run scripts/00 (or supply real data) and scripts/03, scripts/04 first. Missing: ",
                                  paste(need[!file.exists(need)], collapse = ", "))
bm   <- attach.big.matrix(P("features_k3.desc"))          # memory-mapped, ~0 MB in RAM
meta <- fread(P(CFG$meta_file), select = unname(CFG$cols[c("accession", "date", "country", "lineage", "platform")]))
setnames(meta, c("accession", "date", "country", "lineage", "platform"))
meta[, `:=`(month = to_month(date), country = to_country(country))]
meta[, n_frac := bm[, 2]]                                  # row order = feature_index order (same as FASTA)
pca  <- readRDS("results/pca_all.rds")
PLAT <- sort(unique(meta$platform)); LIN <- sort(unique(meta$lineage)); CTRY <- sort(unique(meta$country))
rng  <- range(meta$month, na.rm = TRUE)

ui <- page_sidebar(
  title = "SARS-CoV-2 large dataset: R control panel",
  sidebar = sidebar(width = 300,
    sliderInput("frac", "Sample fraction (%)", 1, 100, 10, step = 1),
    selectInput("lin", "Lineage", LIN, LIN, multiple = TRUE),
    selectInput("ctry", "Country", CTRY, CTRY, multiple = TRUE),
    dateRangeInput("dates", "Collection month", rng[1], rng[2], rng[1], rng[2]),
    radioButtons("backend", "Metadata load backend", c("data.table", "base R", "ff (on disk)")),
    sliderInput("chunk", "Model stream chunk size (rows)", 1000, 20000, 5000, step = 1000),
    radioButtons("mode", "PCA plot mode", c("All points", "Sample", "Hexbin"), "All points"),
    helpText("Synthetic demo data unless SARS2_DATA_DIR points to real data.")),
  navset_card_tab(
    nav_panel("Lineages over time", plotOutput("p_lin", height = 420), textOutput("n_lin")),
    nav_panel("Sequencing quality model",
      actionButton("refit", "Refit (streamed)"), textOutput("fit_note"),
      DTOutput("fit_tbl"), plotOutput("p_q", height = 320)),
    nav_panel("PCA and overplotting", plotOutput("p_pca", height = 460), textOutput("pca_note")),
    nav_panel("Backends and memory",
      layout_columns(actionButton("run_meta", "Time metadata load (selected backend)"),
                     actionButton("run_feat", "Time feature matrix: RAM vs big.matrix")),
      plotOutput("p_mem", height = 240), DTOutput("log_tbl")))
)

server <- function(input, output, session) {
  LOG <- reactiveVal(data.table(step = integer(), action = character(), seconds = numeric(),
                                heap_change_mb = numeric(), heap_total_mb = numeric()))
  log_add <- function(row) isolate(LOG(rbind(LOG(), cbind(step = nrow(LOG()) + 1L,
                       row[, .(action = label, seconds, heap_change_mb = mem_change_mb)],
                       heap_total_mb = round(mem_used() / 1e6, 1)))))

  dat <- reactive({
    req(input$lin, input$ctry)
    m <- measure("filter + sample_frac", {
      x <- meta[lineage %in% input$lin & country %in% input$ctry &
                month >= input$dates[1] & month <= input$dates[2]]
      set.seed(CFG$seed); as.data.table(sample_frac(as.data.frame(x), input$frac / 100))
    })
    log_add(m$row); m$value
  })

  output$p_lin <- renderPlot({
    d <- dat()[!is.na(month), .N, by = .(month, lineage)]
    d <- merge(CJ(month = unique(d$month), lineage = unique(d$lineage)), d, all.x = TRUE)  # zero-fill gaps
    d[is.na(N), N := 0L]; d[, p := N / sum(N), by = month]
    ggplot(d, aes(month, p, fill = lineage)) + geom_area() + theme_minimal(base_size = 13) +
      labs(y = "Share of sampled genomes", x = NULL)
  })
  output$n_lin <- renderText(sprintf("%d genomes after filters and %d%% sampling (full table: %d).",
                                     nrow(dat()), input$frac, nrow(meta)))

  fit <- eventReactive(list(input$refit, input$chunk), {
    prep <- function(d) { d$platform <- factor(d$platform, levels = PLAT); as.data.frame(d) }
    m <- measure(sprintf("stream_fit (chunk %d)", input$chunk),
                 stream_fit(P("model_table.tsv"), n_frac ~ t_years + platform, prep, input$chunk))
    log_add(m$row); m
  }, ignoreNULL = FALSE)
  output$fit_tbl <- renderDT(datatable(fit()$value[, .(term, estimate = signif(estimate, 4), se = signif(se, 3), n)],
                                       options = list(dom = "t"), rownames = FALSE))
  output$fit_note <- renderText(sprintf("engine: %s | %.3f s", fit()$value$engine[1], fit()$row$seconds))
  output$p_q <- renderPlot({
    d <- dat()[!is.na(month), .(n_frac = mean(n_frac)), by = .(month, platform)]
    ggplot(d, aes(month, n_frac, colour = platform)) + geom_line() + theme_minimal(base_size = 13) +
      labs(y = "Mean N fraction", x = NULL)
  })

  output$p_pca <- renderPlot({
    p <- pca[lineage %in% input$lin]
    if (input$mode == "Sample") p <- p[sample(.N, ceiling(.N * input$frac / 100))]
    g <- ggplot(p, aes(PC1, PC2)) + theme_minimal(base_size = 13)
    g <- if (input$mode == "Hexbin") g + geom_hex(bins = 40) else g + geom_point(aes(colour = lineage), size = 0.6, alpha = 0.5)
    t0 <- proc.time()[["elapsed"]]; print(g)
    log_add(data.table(label = sprintf("draw PCA: %s (n=%d)", input$mode, nrow(p)),
                       seconds = round(proc.time()[["elapsed"]] - t0, 3), mem_change_mb = NA_real_))
  })
  output$pca_note <- renderText("Full scatter hides density; hexbin and sampling stay readable as n grows (slides 13-17).")

  observeEvent(input$run_meta, {
    f <- P(CFG$meta_file)
    m <- switch(input$backend,
      "data.table"   = measure("load: data.table fread", fread(f)),
      "base R"       = measure("load: base read.delim", read.delim(f)),
      "ff (on disk)" = measure("load: ff read.table.ffdf", ff::read.table.ffdf(file = f, header = TRUE, sep = "\t", quote = "")))
    log_add(m$row)
  })
  observeEvent(input$run_feat, {
    log_add(measure("features: big.matrix attach (memory-mapped)", attach.big.matrix(P("features_k3.desc")))$row)
    log_add(measure("features: copy all into RAM", { x <- bm[, ]; dim(x) })$row)
  })
  output$log_tbl <- renderDT(datatable(LOG(), options = list(dom = "tp", pageLength = 8), rownames = FALSE))
  output$p_mem <- renderPlot({
    l <- LOG(); req(nrow(l) > 0)
    ggplot(l, aes(step, heap_total_mb)) + geom_line() + geom_point() + theme_minimal(base_size = 13) +
      labs(x = "Action #", y = "R heap (MB)")
  })
}
shinyApp(ui, server)
