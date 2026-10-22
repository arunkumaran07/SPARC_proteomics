base_cfg <- function() {
  cfg <- default_config()
  cfg$input$matrix <- "x.tsv"; cfg$input$sample_sheet <- "s.csv"
  cfg
}

test_that("defaults validate once input files are named", {
  expect_silent(validate_config(base_cfg()))
  expect_error(validate_config(default_config()), "input.matrix")
})

test_that("invalid options are rejected with the allowed choices", {
  cfg <- base_cfg(); cfg$preprocessing$imputation$method <- "magic"
  expect_error(validate_config(cfg), "none, mixed, knn, downshift, mindet")
  cfg <- base_cfg(); cfg$statistics$method <- "anova"
  expect_error(validate_config(cfg), "limma, welch, student")
  cfg <- base_cfg(); cfg$ai$enabled <- "maybe"
  expect_error(validate_config(cfg), "ai.enabled")
})

test_that("species keys and aliases resolve", {
  expect_equal(species_info("Mus musculus")$orgdb, "org.Mm.eg.db")
  expect_equal(species_info("hs")$reactome, "human")
  expect_true(is.na(species_info("other")$reactome))
  expect_error(species_info("unicorn"), "Unknown species")
  cfg <- base_cfg(); cfg$study$species <- "Rattus norvegicus"
  expect_equal(validate_config(cfg)$study$species, "rat")
})

test_that("read_config overrides defaults and resolves relative paths", {
  d <- withr_tempdir()
  writeLines(c("study: {species: mouse}", "input: {matrix: m.tsv, sample_sheet: s.csv}",
               "statistics: {method: welch}"), file.path(d, "c.yml"))
  cfg <- read_config(file.path(d, "c.yml"))
  expect_equal(cfg$statistics$method, "welch")
  expect_equal(cfg$preprocessing$imputation$method, "mixed")  # default kept
  expect_equal(normalizePath(dirname(cfg$input$matrix)), normalizePath(d))
})
