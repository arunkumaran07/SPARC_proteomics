# Benchmark sparcproteomics on a mixed-species dataset with known ratios.
#
# Data: PXD028735 (Van Puyvelde et al., Sci Data 9:126, 2022), Q Exactive HF-X
# DDA, condition A vs B, 3 replicates each, quantified at precursor level by
# WOMBAT-P and distributed with ProteoBench (see benchmark/README.md).
# Expected log2(A/B): human 0, yeast +1, E. coli -2.
#
# Usage (from the repository root, with the package installed):
#   Rscript benchmark/run_benchmark.R          # run and write results
#   Rscript benchmark/run_benchmark.R --check  # also fail if the default
#                                              # settings fall below thresholds

suppressPackageStartupMessages(library(sparcproteomics))
check <- "--check" %in% commandArgs(trailingOnly = TRUE)
expected <- c(HUMAN = 0, YEAST = 1, ECOLI = -2)
data_file <- "benchmark/data/pxd028735_qehfx_dda_precursors.tsv.gz"
stopifnot(file.exists(data_file))
out_dir <- "benchmark/results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# unpack once; the reader takes a plain TSV
tsv <- file.path(tempdir(), "bench_precursors.tsv")
writeLines(readLines(gzfile(data_file)), tsv)
species_of <- local({
  d <- utils::read.delim(tsv, colClasses = c(protein_group = "character", species = "character"))
  tapply(d$species, d$protein_group, `[`, 1)
})
# stable reference for "median_ref": the human proteins (expected ratio 1:1)
ref_file <- file.path(tempdir(), "human_reference.txt")
writeLines(names(species_of)[species_of == "HUMAN"], ref_file)

run_one <- function(summarization, normalization, stat) {
  cfg <- default_config()
  cfg$study$title <- "PXD028735 benchmark"
  cfg$study$species <- "other"
  cfg$input[c("matrix", "format", "sample_sheet", "protein_column", "precursor_column",
              "summarization")] <- list(tsv, "peptide_generic", "benchmark/data/samples.csv",
                                        "protein_group", "precursor", summarization)
  cfg$design$group_column <- "condition"; cfg$design$control <- "B"
  cfg$preprocessing$normalization <- if (normalization == "median_ref") "median" else normalization
  if (normalization == "median_ref") cfg$preprocessing$normalization_reference <- ref_file
  cfg$preprocessing$imputation$method <- "none"   # tests never use imputed values
  cfg$statistics$method <- stat
  cfg$enrichment$enabled <- FALSE; cfg$ai$enabled <- FALSE
  cfg$output$dir <- file.path(tempdir(), paste(summarization, normalization, stat, sep = "_"))
  res <- run_pipeline(cfg, verbose = FALSE)
  t <- res$de$tables[["A_vs_B"]]
  t$species <- unname(species_of[t$protein_id])
  t$expected <- unname(expected[t$species])
  t
}

score <- function(t, adj_p = 0.05, lfc = 0.585) {
  tested <- t[!is.na(t$P.Value), ]
  changed <- tested$species != "HUMAN"
  sig <- !is.na(tested$adj.P.Val) & tested$adj.P.Val < adj_p & abs(tested$logFC) >= lfc
  correct_dir <- sign(tested$logFC) == sign(tested$expected)
  tp <- sig & changed & correct_dir
  fp <- sig & !tp
  sig_p <- !is.na(tested$adj.P.Val) & tested$adj.P.Val < adj_p
  med <- tapply(tested$logFC, tested$species, stats::median)
  mad <- tapply(tested$logFC, tested$species, stats::mad)
  data.frame(
    tested_human = sum(tested$species == "HUMAN"),
    tested_yeast = sum(tested$species == "YEAST"),
    tested_ecoli = sum(tested$species == "ECOLI"),
    significant = sum(sig),
    sensitivity = sum(tp) / sum(changed),
    fdp = if (sum(sig)) sum(fp) / sum(sig) else 0,
    fdp_adjp_only = if (sum(sig_p)) sum(sig_p & !(changed & correct_dir)) / sum(sig_p) else 0,
    bias_human = med[["HUMAN"]] - 0, bias_yeast = med[["YEAST"]] - 1, bias_ecoli = med[["ECOLI"]] + 2,
    mad_human = mad[["HUMAN"]], mad_yeast = mad[["YEAST"]], mad_ecoli = mad[["ECOLI"]])
}

grid <- expand.grid(stat = c("limma", "welch", "student"),
                    normalization = c("median", "median_ref", "none"),
                    summarization = c("median_polish", "top3", "sum"),
                    stringsAsFactors = FALSE)
results <- list(); tables <- list()
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  t <- run_one(g$summarization, g$normalization, g$stat)
  key <- paste(g$summarization, g$normalization, g$stat, sep = " / ")
  tables[[key]] <- t
  results[[key]] <- cbind(g[c("summarization", "normalization", "stat")], score(t))
  message(sprintf("%-32s sensitivity %.3f  FDP %.3f", key, results[[key]]$sensitivity,
                  results[[key]]$fdp))
}
metrics <- do.call(rbind, results); rownames(metrics) <- NULL
utils::write.csv(metrics, file.path(out_dir, "metrics.csv"), row.names = FALSE)

# figure for the default settings
default_key <- "median_polish / median / limma"
if (requireNamespace("ggplot2", quietly = TRUE)) {
  d <- tables[[default_key]]
  d <- d[!is.na(d$logFC), ]
  d$species <- factor(d$species, levels = c("HUMAN", "YEAST", "ECOLI"))
  cols <- c(HUMAN = "#7F7F7F", YEAST = "#E69F00", ECOLI = "#0072B2")
  p1 <- ggplot2::ggplot(d, ggplot2::aes(AveExpr, logFC, colour = species)) +
    ggplot2::geom_point(size = 0.6, alpha = 0.5) +
    ggplot2::geom_hline(yintercept = expected, colour = cols[names(expected)], linetype = 2) +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::labs(x = "Mean log2 intensity", y = "log2(A/B)", colour = NULL,
                  title = "Measured vs expected ratios (dashed)") +
    ggplot2::theme_bw(base_size = 10)
  p2 <- ggplot2::ggplot(d, ggplot2::aes(species, logFC, fill = species)) +
    ggplot2::geom_boxplot(outlier.size = 0.3, show.legend = FALSE) +
    ggplot2::geom_point(data = data.frame(species = factor(names(expected), levels(d$species)),
                                          logFC = expected),
                        shape = 23, size = 2.5, fill = "white", show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::labs(x = NULL, y = "log2(A/B)", title = "Distribution per species (diamond = expected)") +
    ggplot2::theme_bw(base_size = 10)
  fig <- file.path(out_dir, "default_settings.png")
  grDevices::png(fig, width = 2400, height = 1000, res = 220)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2, widths = c(1.4, 1))))
  print(p1, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p2, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grDevices::dev.off()
}

fmt <- function(x, d = 3) formatC(as.numeric(x), format = "f", digits = d)
pct <- function(x) paste0(formatC(100 * as.numeric(x), format = "f", digits = 1), "%")
row_md <- function(r, cols) paste0("| ", paste(vapply(cols, function(k) switch(k,
  sensitivity = pct(r[[k]]), fdp = pct(r[[k]]), fdp_adjp_only = pct(r[[k]]),
  bias_human = fmt(r[[k]], 2), bias_yeast = fmt(r[[k]], 2), bias_ecoli = fmt(r[[k]], 2),
  as.character(r[[k]])), ""), collapse = " | "), " |")
tab_md <- function(df, cols, header) c(
  paste0("| ", paste(header, collapse = " | "), " |"),
  paste0("|", paste(rep("---", length(header)), collapse = "|"), "|"),
  vapply(seq_len(nrow(df)), function(i) row_md(df[i, ], cols), ""))
nm <- c(median = "median (all proteins)", median_ref = "median (human reference)", none = "none")
metrics$normalization_label <- unname(nm[metrics$normalization])
mp <- metrics[metrics$summarization == "median_polish", ]
lm <- metrics[metrics$stat == "limma", ]
dflt <- metrics[metrics$summarization == "median_polish" & metrics$normalization == "median" &
                  metrics$stat == "limma", ]
md <- c(
  "# Benchmark results",
  "",
  sprintf("Generated %s with sparcproteomics %s. Rerun with `Rscript benchmark/run_benchmark.R`.",
          format(Sys.Date()), utils::packageVersion("sparcproteomics")),
  "",
  "**Dataset:** PXD028735 (Van Puyvelde et al. 2022), Q Exactive HF-X DDA, condition A vs B,",
  "3 replicates each; precursor-level quantification from ProteoBench (WOMBAT-P).",
  sprintf("%d yeast, %d *E. coli* and %d human proteins tested with default settings.",
          dflt$tested_yeast, dflt$tested_ecoli, dflt$tested_human),
  "Expected log2(A/B): human 0, yeast +1, *E. coli* -2.",
  "",
  "**Scoring:** significant = BH-adjusted p < 0.05 and |log2FC| >= 0.585. A true positive is a",
  "yeast or *E. coli* protein called significant in the expected direction; any other significant",
  "call is false. Sensitivity = true positives / yeast + *E. coli* proteins tested. FDP = false /",
  "all significant calls. \"FDP, adj. p only\" drops the fold-change cutoff. Bias = median log2FC",
  "minus the expected value. No imputation is used for testing.",
  "",
  "## Effect of the statistical test (median-polish summarization)",
  "",
  tab_md(mp, c("normalization_label", "stat", "significant", "sensitivity", "fdp", "fdp_adjp_only"),
         c("Normalization", "Test", "Significant", "Sensitivity", "FDP", "FDP, adj. p only")),
  "",
  "## Effect of summarization and normalization (limma)",
  "",
  tab_md(lm, c("summarization", "normalization_label", "sensitivity", "fdp", "bias_human",
               "bias_yeast", "bias_ecoli"),
         c("Summarization", "Normalization", "Sensitivity", "FDP", "Bias human",
           "Bias yeast", "Bias E. coli")),
  "",
  "All 27 combinations, with precision (MAD) per species, are in `results/metrics.csv`.",
  "",
  "![Default settings](results/default_settings.png)")
writeLines(md, "benchmark/RESULTS.md")
message("Wrote benchmark/RESULTS.md and ", out_dir)

utils::write.csv(tables[[default_key]], gzfile(file.path(out_dir, "DE_default_settings.csv.gz")),
                 row.names = FALSE)

# Regression guards, set just below the results of sparcproteomics 0.1.0 so
# that a code change that makes the default analysis worse fails CI. They are
# not ideal targets: the E. coli ratio is compressed in every setting,
# including no normalization, so that bias comes from the data.
if (check) {
  m <- metrics[metrics$summarization == "median_polish" & metrics$normalization == "median" &
                 metrics$stat == "limma", ]
  r <- metrics[metrics$summarization == "median_polish" & metrics$normalization == "median_ref" &
                 metrics$stat == "limma", ]
  fails <- c(
    if (m$sensitivity < 0.75) sprintf("default sensitivity %.3f < 0.75", m$sensitivity),
    if (m$fdp > 0.02) sprintf("default FDP %.3f > 0.02", m$fdp),
    if (abs(m$bias_human) > 0.20) sprintf("default |human bias| %.2f > 0.20", m$bias_human),
    if (abs(m$bias_yeast) > 0.35) sprintf("default |yeast bias| %.2f > 0.35", m$bias_yeast),
    if (abs(m$bias_ecoli) > 0.75) sprintf("default |E. coli bias| %.2f > 0.75", m$bias_ecoli),
    if (r$sensitivity < 0.83) sprintf("reference-normalized sensitivity %.3f < 0.83", r$sensitivity),
    if (abs(r$bias_human) > 0.05) sprintf("reference-normalized |human bias| %.2f > 0.05", r$bias_human))
  if (length(fails)) stop("Benchmark check failed: ", paste(fails, collapse = "; "), call. = FALSE)
  message("Benchmark check passed: default sensitivity ", fmt(m$sensitivity), ", FDP ", fmt(m$fdp))
}
