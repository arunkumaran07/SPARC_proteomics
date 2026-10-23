test_that("median polish recovers sample effects from an additive model", {
  set.seed(1)
  truth <- c(20, 21, 22.5, 19)                  # protein abundance per sample
  offset <- c(0, -3, 2, -1.5, 4)                # precursor response factors
  x <- outer(offset, truth, "+") + matrix(rnorm(20, 0, 0.01), 5)
  x[2, 3] <- NA; x[5, 1] <- NA
  x[4, 2] <- x[4, 2] + 6                        # one outlying value
  colnames(x) <- paste0("s", 1:4)
  est <- summarize_precursors(x, rep("P1", 5), "median_polish")$x[1, ]
  expect_equal(unname(est - est[1]), truth - truth[1], tolerance = 0.05)
})

test_that("top3 and sum behave as documented", {
  x <- log2(rbind(c(100, 200), c(10, 20), c(1000, 2000), c(1, NA)))
  colnames(x) <- c("a", "b")
  top3 <- summarize_precursors(x, rep("P", 4), "top3")$x
  expect_equal(unname(top3[1, ]), log2(c(mean(c(1000, 100, 10)), mean(c(2000, 200, 20)))))
  s <- summarize_precursors(x, rep("P", 4), "sum")$x
  expect_equal(unname(s[1, ]), log2(c(1111, 2220)))
})

test_that("proteins unseen in a sample stay missing and min_precursors filters", {
  x <- rbind(c(20, NA, 21), c(22, NA, 23), c(18, 19, 20))
  colnames(x) <- c("a", "b", "c")
  r <- summarize_precursors(x, c("P1", "P1", "P2"))
  expect_true(is.na(r$x["P1", "b"]))
  expect_equal(unname(r$n_precursors), c(2L, 1L))
  expect_equal(rownames(summarize_precursors(x, c("P1", "P1", "P2"), min_precursors = 2)$x), "P1")
})

test_that("DIA-NN main report is read, q-filtered and pivoted", {
  f <- test_path("fixtures", "diann_report_subset.tsv")
  raw <- utils::read.delim(f)
  tab <- read_precursor_table(f, "diann_report")
  ok <- raw$Q.Value <= 0.01 & raw$PG.Q.Value <= 0.01
  expect_equal(ncol(tab$intensities), length(unique(raw$Run)))
  expect_equal(nrow(tab$intensities), length(unique(raw$Precursor.Id[ok])))
  # a value lands in the right cell
  i <- which(ok)[10]
  expect_equal(tab$intensities[raw$Precursor.Id[i], raw$Run[i]], raw$Precursor.Normalised[i])
  expect_equal(tab$features$protein_id[tab$features$precursor_id == raw$Precursor.Id[i]],
               raw$Protein.Group[i])
})

test_that("pipeline runs from a DIA-NN main report", {
  d <- withr_tempdir()
  f <- test_path("fixtures", "diann_report_subset.tsv")
  runs <- unique(utils::read.delim(f)$Run)
  utils::write.csv(data.frame(sample_id = runs,
                              condition = ifelse(grepl("Condition_A", runs), "A", "B")),
                   file.path(d, "samples.csv"), row.names = FALSE)
  cfg <- default_config()
  cfg$input$matrix <- f; cfg$input$format <- "diann_report"
  cfg$input$sample_sheet <- file.path(d, "samples.csv")
  cfg$design$group_column <- "condition"; cfg$design$control <- "B"
  cfg$output$dir <- file.path(d, "results")
  cfg$enrichment$enabled <- FALSE; cfg$ai$enabled <- FALSE
  res <- run_pipeline(cfg, verbose = FALSE)
  expect_equal(ncol(res$se), 6)
  expect_true("n_precursors" %in% names(SummarizedExperiment::rowData(res$se)))
  expect_true(file.exists(file.path(d, "results", "DE_A_vs_B.csv")))
})

test_that("peptide_generic input with simulated precursors recovers changes", {
  set.seed(5)
  d <- withr_tempdir()
  n_prot <- 400; per <- 4
  prot <- rep(sprintf("P%03d", seq_len(n_prot)), each = per)
  ab <- rep(rnorm(n_prot, 22, 2), each = per) + rnorm(n_prot * per, 0, 1.5)
  x <- ab + matrix(rnorm(n_prot * per * 8, 0, 0.3), ncol = 8)
  changed <- sprintf("P%03d", 1:40)
  x[prot %in% changed, 5:8] <- x[prot %in% changed, 5:8] + 1.5
  x[matrix(runif(length(x)) < 0.05, nrow(x))] <- NA
  df <- data.frame(peptide = paste0("pep", seq_along(prot)), protein = prot, 2^x)
  names(df)[3:10] <- paste0("S", 1:8)
  utils::write.table(df, file.path(d, "peptides.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.csv(data.frame(sample_id = paste0("S", 1:8), group = rep(c("Ctrl", "Trt"), each = 4)),
                   file.path(d, "samples.csv"), row.names = FALSE)
  for (m in c("median_polish", "top3", "sum")) {
    cfg <- default_config()
    cfg$input[c("matrix", "format", "sample_sheet", "protein_column", "precursor_column",
                "summarization")] <- list(file.path(d, "peptides.tsv"), "peptide_generic",
                                          file.path(d, "samples.csv"), "protein", "peptide", m)
    cfg$design$control <- "Ctrl"; cfg$output$dir <- file.path(d, m)
    cfg$enrichment$enabled <- FALSE; cfg$ai$enabled <- FALSE
    res <- run_pipeline(cfg, verbose = FALSE)
    t <- res$de$tables[[1]]
    sig <- t$protein_id[t$significant]
    expect_gt(mean(changed %in% sig), 0.8)
    expect_lt(mean(!sig %in% changed), 0.1)
  }
})

test_that("peptide_generic requires a protein column", {
  cfg <- default_config()
  cfg$input$matrix <- "x"; cfg$input$sample_sheet <- "y"; cfg$input$format <- "peptide_generic"
  expect_error(validate_config(cfg), "protein_column")
})

test_that("reference-set median normalization removes the shift when many proteins change", {
  set.seed(9)
  n <- 200
  x <- matrix(rnorm(n * 6, 22, 0.1), n) + rnorm(n, 0, 2)
  x[1:90, 4:6] <- x[1:90, 4:6] + 2               # 45% of proteins rise in group B
  rownames(x) <- paste0("P", seq_len(n)); colnames(x) <- paste0("S", 1:6)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(log2 = x),
    rowData = S4Vectors::DataFrame(protein_id = rownames(x), gene = NA_character_),
    colData = S4Vectors::DataFrame(group = factor(rep(c("A", "B"), each = 3))))
  stable <- rownames(x)[91:n]
  fc <- function(s) { y <- SummarizedExperiment::assay(s, "norm")[stable, ]
    stats::median(rowMeans(y[, 4:6]) - rowMeans(y[, 1:3])) }
  expect_gt(abs(fc(pp_normalize(se, "median"))), 0.5)          # biased by the changers
  expect_lt(abs(fc(pp_normalize(se, "median", reference = stable))), 0.05)
  expect_error(pp_normalize(se, "median", reference = c("P1", "P2")), "at least 20")
})
