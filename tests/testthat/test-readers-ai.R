test_that("annotated workbook is read, including duplicate log2 columns", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")
  set.seed(3)
  ids <- rep(sprintf("M%d", 1:6), each = 2)            # linear + log2 per sample
  grp <- rep(rep(c("Ctrl", "Drug"), each = 3), each = 2)
  lin <- matrix(2^rnorm(50 * 12, 20, 1), 50)
  lin[, seq(2, 12, 2)] <- log2(lin[, seq(1, 11, 2)])   # second copy is log2
  top <- rbind(c("Protein", ids), c("Group", grp), c("Sex", rep(c("F", "M"), 6)))
  body <- cbind(sprintf("G%02d", 1:50), format(lin, digits = 8))
  f <- tempfile(fileext = ".xlsx")
  openxlsx::write.xlsx(as.data.frame(rbind(top, body)), f, colNames = FALSE)

  tab <- read_annotated_xlsx(f)
  expect_equal(names(tab$samples), c("sample_id", "Group", "Sex"))
  expect_equal(ncol(tab$intensities), 12)

  cfg <- default_config()
  cfg$input$matrix <- f; cfg$input$format <- "annotated_xlsx"
  cfg$input$duplicate_columns <- "log2"
  cfg$design$group_column <- "Group"; cfg$design$control <- "Ctrl"
  se <- build_dataset(cfg)
  expect_equal(ncol(se), 6)
  expect_true(stats::median(SummarizedExperiment::assay(se, "log2"), na.rm = TRUE) < 25)

  cfg$input$duplicate_columns <- "error"
  expect_error(build_dataset(cfg), "Repeated sample names")
})

test_that("AI is used only when enabled and a key is present", {
  cfg <- default_config()
  cfg$ai$key_env <- "SP_TEST_NO_SUCH_KEY"
  expect_false(ai_status(cfg)$use)
  cfg$ai$enabled <- FALSE
  expect_equal(ai_status(cfg)$reason, "disabled in config")
  skip_if_not_installed("httr2")
  Sys.setenv(SP_TEST_KEY = "dummy"); on.exit(Sys.unsetenv("SP_TEST_KEY"))
  cfg$ai$key_env <- "SP_TEST_KEY"; cfg$ai$enabled <- "auto"
  expect_true(ai_status(cfg)$use)
})

test_that("AI prompt carries a null-result guard and no raw data", {
  d <- withr_tempdir()
  cfg <- read_config(simulate_dataset(d, frac_de = 0, n_proteins = 300))
  res <- run_pipeline(cfg, verbose = FALSE)
  p <- sparcproteomics:::.ai_prompts(res$de, res$se, cfg, NULL, 10)
  expect_match(p$user, "null result")
  expect_false(grepl("S01", p$user))   # sample-level values are never sent
})
