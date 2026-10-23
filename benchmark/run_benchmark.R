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

# the package defaults this benchmark is scored against
dc <- default_config()
DEFAULT <- list(summarization = dc$input$summarization,
                normalization = dc$preprocessing$normalization,
                stat = if (dc$statistics$significance == "treat") "limma_treat" else dc$statistics$method)
key_of <- function(s, n, t) paste(s, n, t, sep = " / ")
default_key <- key_of(DEFAULT$summarization, DEFAULT$normalization, DEFAULT$stat)

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
  cfg$statistics$method <- if (stat == "limma_treat") "limma" else stat
  cfg$statistics$significance <- if (stat == "limma_treat") "treat" else "cutoff"
  cfg$statistics$treat_fc <- 1.5
  cfg$enrichment$enabled <- FALSE; cfg$ai$enabled <- FALSE
  cfg$output$dir <- file.path(tempdir(), paste(summarization, normalization, stat, sep = "_"))
  res <- run_pipeline(cfg, verbose = FALSE)
  t <- res$de$tables[["A_vs_B"]]
  t$species <- unname(species_of[t$protein_id])
  t$expected <- unname(expected[t$species])
  det <- res$detection$tables[["A_vs_B"]]
  det$species <- unname(species_of[det$protein_id])
  list(de = t, detection = det)
}

# Uses the pipeline's own significance call, so each rule (cutoff or TREAT)
# is scored as it would be reported.
score <- function(t, adj_p = 0.05) {
  tested <- t[!is.na(t$P.Value), ]
  changed <- tested$species != "HUMAN"
  sig <- tested$significant
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

# On/off proteins: correct when a yeast protein is only in A (yeast is higher
# in A) or an E. coli protein only in B; human on/off calls are false.
score_detection <- function(d) {
  on <- nzchar(d$on_off)
  strict <- on & ((d$on_off == "only_in_A" & d$fraction_A == 1) |
                  (d$on_off == "only_in_B" & d$fraction_B == 1))
  ok <- (d$on_off == "only_in_A" & d$species == "YEAST") |
        (d$on_off == "only_in_B" & d$species == "ECOLI")
  data.frame(rule = c(">= 2 of 3 in one group, 0 in the other", "3 of 3 in one group, 0 in the other"),
             listed = c(sum(on), sum(strict)),
             correct = c(sum(on & ok), sum(strict & ok)),
             human = c(sum(on & d$species == "HUMAN"), sum(strict & d$species == "HUMAN")),
             wrong_direction = c(sum(on & !ok & d$species != "HUMAN"),
                                 sum(strict & !ok & d$species != "HUMAN")))
}

grid <- expand.grid(stat = c("limma", "limma_treat", "welch", "student"),
                    normalization = c("median", "mode", "median_ref", "none"),
                    summarization = c("median_polish", "maxlfq", "top3", "sum"),
                    stringsAsFactors = FALSE)
results <- list(); tables <- list(); detections <- list()
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  key <- key_of(g$summarization, g$normalization, g$stat)
  r <- run_one(g$summarization, g$normalization, g$stat)
  tables[[key]] <- r$de
  if (g$stat == "limma" && g$normalization == "median") detections[[g$summarization]] <- r$detection
  results[[key]] <- cbind(g[c("summarization", "normalization", "stat")], score(r$de))
  message(sprintf("%-40s sensitivity %.3f  FDP %.3f", key, results[[key]]$sensitivity,
                  results[[key]]$fdp))
}
metrics <- do.call(rbind, results); rownames(metrics) <- NULL
utils::write.csv(metrics, file.path(out_dir, "metrics.csv"), row.names = FALSE)
det_score <- score_detection(detections[[DEFAULT$summarization]])
utils::write.csv(det_score, file.path(out_dir, "detection.csv"), row.names = FALSE)

# figure for the default settings
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
nm <- c(median = "median (all proteins)", mode = "mode of ratios",
        median_ref = "median (human reference)", none = "none")
st <- c(limma = "limma, p + 1.5-fold cutoff", limma_treat = "limma TREAT, 1.5-fold",
        welch = "Welch, p + 1.5-fold cutoff", student = "Student, p + 1.5-fold cutoff")
metrics$normalization_label <- unname(nm[metrics$normalization])
metrics$stat_label <- unname(st[metrics$stat])
dflt <- metrics[metrics$summarization == DEFAULT$summarization &
                  metrics$normalization == DEFAULT$normalization & metrics$stat == DEFAULT$stat, ]
t1 <- metrics[metrics$summarization == DEFAULT$summarization &
                metrics$normalization %in% c("median", "mode"), ]
t1 <- t1[order(match(t1$normalization, c("median", "mode")), match(t1$stat, names(st))), ]
t2 <- metrics[metrics$stat == "limma", ]
t2 <- t2[order(match(t2$summarization, c("median_polish", "maxlfq", "top3", "sum")),
               match(t2$normalization, names(nm))), ]
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
  sprintf("**Defaults scored:** %s summarization, %s normalization, %s.",
          DEFAULT$summarization, nm[[DEFAULT$normalization]], st[[DEFAULT$stat]]),
  "",
  "**Scoring:** significant as reported by the pipeline: BH-adjusted p < 0.05 and |log2FC| >= 0.585",
  "(1.5-fold), or for TREAT, BH-adjusted p < 0.05 for the hypothesis |FC| > 1.5. A true positive is",
  "a yeast or *E. coli* protein called significant in the expected direction; any other",
  "significant call is false. Sensitivity = true positives / yeast + *E. coli* proteins tested.",
  "FDP = false / all significant calls. \"FDP, adj. p only\" ignores the fold-change requirement.",
  "Bias = median log2FC minus the expected value. No imputation is used for testing.",
  "",
  sprintf("## Significance rule and test (%s summarization)", DEFAULT$summarization),
  "",
  tab_md(t1, c("normalization_label", "stat_label", "significant", "sensitivity", "fdp", "fdp_adjp_only"),
         c("Normalization", "Test and rule", "Significant", "Sensitivity", "FDP", "FDP, adj. p only")),
  "",
  "## Summarization and normalization (limma, p + 1.5-fold cutoff)",
  "",
  tab_md(t2, c("summarization", "normalization_label", "sensitivity", "fdp", "fdp_adjp_only",
               "bias_human", "bias_yeast", "bias_ecoli"),
         c("Summarization", "Normalization", "Sensitivity", "FDP", "FDP, adj. p only",
           "Bias human", "Bias yeast", "Bias E. coli")),
  "",
  "## Presence/absence (on/off) proteins",
  "",
  "With 3 vs 3 samples Fisher's exact test cannot reach significance (smallest possible p = 0.1),",
  "so on/off proteins are listed descriptively. Correct = yeast present only in A, or *E. coli*",
  "present only in B.",
  "",
  tab_md(det_score, c("rule", "listed", "correct", "human", "wrong_direction"),
         c("Rule", "Listed", "Correct", "Human (false)", "Wrong direction")),
  "",
  if (all(det_score$listed == 0)) paste(
    "**Not evaluable on this dataset:** every protein is quantified in at least one run of each",
    "condition, so there are no on/off proteins to find. The module is covered by unit tests;",
    "benchmarking it needs a dataset with proteins absent from a whole condition.") else "",
  "",
  sprintf("All %d combinations, with precision (MAD) per species, are in `results/metrics.csv`.",
          nrow(metrics)),
  "",
  "![Default settings](results/default_settings.png)")
writeLines(md, "benchmark/RESULTS.md")
message("Wrote benchmark/RESULTS.md and ", out_dir)

utils::write.csv(tables[[default_key]], gzfile(file.path(out_dir, "DE_default_settings.csv.gz")),
                 row.names = FALSE)

# Regression guards for the default settings, set just below the results of
# the version that introduced them, so a change that makes the default
# analysis worse fails CI. They are not ideal targets: the E. coli ratio is
# compressed in every setting, including no normalization (a property of the
# data).
THRESHOLDS <- list(sensitivity = 0.80, fdp = 0.02, bias_human = 0.10, bias_yeast = 0.25,
                   bias_ecoli = 0.70)
if (check) {
  m <- dflt
  fails <- c(
    if (m$sensitivity < THRESHOLDS$sensitivity)
      sprintf("default sensitivity %.3f < %.2f", m$sensitivity, THRESHOLDS$sensitivity),
    if (m$fdp > THRESHOLDS$fdp) sprintf("default FDP %.3f > %.2f", m$fdp, THRESHOLDS$fdp),
    if (abs(m$bias_human) > THRESHOLDS$bias_human)
      sprintf("default |human bias| %.2f > %.2f", m$bias_human, THRESHOLDS$bias_human),
    if (abs(m$bias_yeast) > THRESHOLDS$bias_yeast)
      sprintf("default |yeast bias| %.2f > %.2f", m$bias_yeast, THRESHOLDS$bias_yeast),
    if (abs(m$bias_ecoli) > THRESHOLDS$bias_ecoli)
      sprintf("default |E. coli bias| %.2f > %.2f", m$bias_ecoli, THRESHOLDS$bias_ecoli))
  if (length(fails)) stop("Benchmark check failed: ", paste(fails, collapse = "; "), call. = FALSE)
  message("Benchmark check passed: default sensitivity ", fmt(m$sensitivity), ", FDP ", fmt(m$fdp))
}
