run_sim <- function(..., edit = identity, sim = list()) {
  d <- withr_tempdir()
  cp <- do.call(simulate_dataset, c(list(dir = d), sim))
  cfg <- edit(read_config(cp))
  res <- run_pipeline(cfg, verbose = FALSE)
  list(res = res, rec = recovery(res, file.path(d, "truth.csv")), dir = d)
}

test_that("each statistical method recovers simulated changes", {
  for (m in c("limma", "welch", "student")) {
    r <- run_sim(edit = function(c) { c$statistics$method <- m; c },
                 sim = list(n_groups = 3, n_per_group = 5))
    expect_equal(r$res$de$method, m)
    # unmoderated t-tests are expected to have less power at small n
    expect_true(all(r$rec[, "power"] > if (m == "limma") 0.7 else 0.4), info = m)
    expect_true(all(r$rec[, "fdp"] < 0.15), info = m)
  }
})

test_that("every imputation and normalization option runs", {
  for (im in c("none", "mixed", "knn", "downshift", "mindet")) {
    r <- run_sim(edit = function(c) { c$preprocessing$imputation$method <- im; c })
    imp <- SummarizedExperiment::assay(r$res$se, "imputed")
    if (im == "none") expect_true(anyNA(imp)) else expect_false(anyNA(imp), info = im)
    # tests use the un-imputed matrix, so DE results do not depend on imputation
    expect_true(all(r$rec[, "power"] > 0.7), info = im)
  }
  for (nm in c("none", "median", "quantile", "cyclicloess")) {
    r <- run_sim(edit = function(c) { c$preprocessing$normalization <- nm; c })
    expect_equal(S4Vectors::metadata(r$res$se)$normalization, nm)
  }
})

test_that("control choice and custom comparisons are respected", {
  r <- run_sim(sim = list(n_groups = 3), edit = function(c) {
    c$design$control <- "TreatmentA"
    c$design$comparisons <- list("TreatmentB vs TreatmentA", "Control - TreatmentA")
    c
  })
  expect_equal(levels(r$res$se$group)[1], "TreatmentA")
  expect_equal(names(r$res$de$tables), c("TreatmentB_vs_TreatmentA", "Control_vs_TreatmentA"))
})

test_that("limma reports a global F-test for >2 groups and accepts covariates", {
  r <- run_sim(sim = list(n_groups = 4), edit = function(c) {
    c$design$covariates <- list("sex"); c })
  expect_false(is.null(r$res$de$f_test))
  expect_true("sex_M" %in% colnames(r$res$de$design))
  expect_equal(nrow(r$res$de$summary), 3)
})

test_that("proteins with too few values in a group are not tested", {
  r <- run_sim(edit = function(c) { c$statistics$min_valid_per_group <- 4; c })
  t <- r$res$de$tables[[1]]
  low <- t$n_valid_TreatmentA < 4 | t$n_valid_Control < 4
  expect_true(any(low))
  expect_true(all(is.na(t$P.Value[low])))
})

test_that("all simulated formats read identically", {
  n <- sapply(c("diann", "spectronaut", "maxquant"), function(f) {
    r <- run_sim(sim = list(format = f, seed = 7))
    r$res$de$summary$significant
  })
  expect_equal(unname(n["diann"]), unname(n["spectronaut"]))
  expect_true(abs(n["maxquant"] - n["diann"]) <= 3)  # 5 decoy rows are removed
})

test_that("missing control or design columns stop the run", {
  d <- withr_tempdir(); cfg <- read_config(simulate_dataset(d))
  cfg$design$control <- "Placebo"
  expect_error(run_pipeline(cfg, verbose = FALSE), "Control group 'Placebo' not found")
  cfg <- read_config(simulate_dataset(d)); cfg$design$covariates <- list("age")
  expect_error(run_pipeline(cfg, verbose = FALSE), "not in the sample annotations")
})
