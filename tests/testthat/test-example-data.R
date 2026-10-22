test_that("DIA-NN example runs end to end and recovers the planted changes", {
  d <- withr_tempdir()
  cfg <- example_data("diann", d)
  res <- run_pipeline(cfg, verbose = FALSE)

  # all 15 runs matched to the sample sheet despite Windows paths and .d suffix
  expect_equal(ncol(res$se), 15)
  expect_equal(levels(res$se$group), c("Control", "TreatmentA", "TreatmentB"))
  expect_true("batch_B2" %in% colnames(res$de$design))
  expect_false(is.null(res$de$f_test))

  truth <- utils::read.csv(file.path(d, "truth.csv"))
  for (g in c("TreatmentA", "TreatmentB")) {
    t <- res$de$tables[[paste0(g, "_vs_Control")]]
    tp <- truth$protein_id[truth$group == g]
    sig <- t$protein_id[t$significant]
    expect_gt(mean(tp %in% sig), 0.75)
    expect_lt(mean(!sig %in% tp), 0.10)
    # directions agree with what was planted
    m <- merge(t[t$significant, ], truth[truth$group == g, ], by = "protein_id")
    expect_true(all(sign(m$logFC) == sign(m$true_log2fc)))
  }
  # the coherent gene set dominates the TreatmentA hits
  etc <- truth$protein_id[grepl("electron transport", truth$set)]
  top <- utils::head(res$de$tables$TreatmentA_vs_Control$protein_id, 40)
  expect_gt(mean(top %in% etc), 0.4)

  for (f in c("DE_summary.csv", "F_test.csv", "PCA.png", "run_log.txt", "results.rds"))
    expect_true(file.exists(file.path(d, "results", f)), info = f)
})

test_that("annotated-workbook example reads duplicate log2 columns correctly", {
  skip_if_not_installed("readxl")
  d <- withr_tempdir()
  res <- run_pipeline(example_data("annotated", d), verbose = FALSE)
  expect_equal(ncol(res$se), 10)                      # 20 columns -> 10 samples
  expect_equal(levels(res$se$group), c("Control", "TreatmentA"))
  truth <- utils::read.csv(file.path(d, "truth.csv"))
  t <- res$de$tables$TreatmentA_vs_Control
  # the workbook carries gene symbols only, so compare by gene
  expect_gt(mean(truth$gene %in% t$gene[t$significant]), 0.75)
})

test_that("both examples give the same TreatmentA result from different inputs", {
  skip_if_not_installed("readxl")
  d1 <- withr_tempdir(); d2 <- withr_tempdir()
  c1 <- read_config(example_data("diann", d1))
  c1$design$comparisons <- list("TreatmentA vs Control"); c1$design$covariates <- list()
  r1 <- run_pipeline(c1, verbose = FALSE)
  r2 <- run_pipeline(example_data("annotated", d2), verbose = FALSE)
  a <- r1$de$tables[[1]]; b <- r2$de$tables[[1]]
  m <- merge(a, b, by = "gene")
  expect_gt(stats::cor(m$logFC.x, m$logFC.y, use = "complete.obs"), 0.95)
})
