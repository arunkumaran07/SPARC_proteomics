make_se <- function(x, group) {
  SummarizedExperiment::SummarizedExperiment(
    assays = list(log2 = x, norm = x),
    rowData = S4Vectors::DataFrame(protein_id = rownames(x), gene = NA_character_),
    colData = S4Vectors::DataFrame(group = factor(group, levels = unique(group))))
}

test_that("MaxLFQ recovers sample effects and preserves total intensity", {
  set.seed(2)
  truth <- c(20, 21, 19.5, 22, 20.5)
  offset <- c(0, -2, 3, -1)
  x <- outer(offset, truth, "+") + matrix(rnorm(20, 0, 0.01), 4)
  x[1, 2] <- NA; x[3, 5] <- NA; x[2, 1] <- NA
  colnames(x) <- paste0("s", 1:5)
  est <- summarize_precursors(x, rep("P", 4), "maxlfq")$x[1, ]
  expect_equal(unname(est - est[1]), truth - truth[1], tolerance = 0.03)
  expect_equal(sum(2^est), sum(2^x, na.rm = TRUE), tolerance = 1e-6)
})

test_that("MaxLFQ leaves samples unlinked to the main group missing", {
  x <- rbind(c(20, 21, NA, NA), c(18, 19, NA, NA), c(NA, NA, 25, NA))
  colnames(x) <- paste0("s", 1:4)
  est <- summarize_precursors(x, rep("P", 3), "maxlfq")$x[1, ]
  expect_false(anyNA(est[1:2]))
  expect_true(all(is.na(est[3:4])))
})

test_that("mode normalization is unbiased when many proteins change one way", {
  set.seed(3)
  n <- 1000
  x <- matrix(rnorm(n * 6, 0, 0.15), n) + rnorm(n, 22, 2)
  x[, 4:6] <- x[, 4:6] + rep(c(0.3, -0.2, 0.1), each = n)   # loading differences
  x[1:350, 4:6] <- x[1:350, 4:6] + 1.5                      # 35% rise in group B
  rownames(x) <- paste0("P", seq_len(n)); colnames(x) <- paste0("S", 1:6)
  se <- make_se(x, rep(c("A", "B"), each = 3))
  fc_unchanged <- function(s) {
    y <- SummarizedExperiment::assay(s, "norm")[351:n, ]
    stats::median(rowMeans(y[, 4:6]) - rowMeans(y[, 1:3]))
  }
  expect_gt(abs(fc_unchanged(pp_normalize(se, "median"))), 0.3)
  expect_lt(abs(fc_unchanged(pp_normalize(se, "mode"))), 0.05)
  # with balanced changes, mode and median agree
  xb <- x; xb[1:350, 4:6] <- xb[1:350, 4:6] - 1.5
  xb[1:175, 4:6] <- xb[1:175, 4:6] + 1.5; xb[176:350, 4:6] <- xb[176:350, 4:6] - 1.5
  sb <- make_se(xb, rep(c("A", "B"), each = 3))
  expect_lt(abs(fc_unchanged(pp_normalize(sb, "mode")) - fc_unchanged(pp_normalize(sb, "median"))), 0.05)
})

test_that("the centering check warns about an off-centre log2FC distribution", {
  set.seed(8); sparcproteomics:::log_reset()
  de <- list(tables = list(A_vs_B = data.frame(logFC = c(rnorm(500, -0.3, 0.1), rnorm(100, 2, 0.2)))))
  sparcproteomics:::.check_centering(de, "median")
  expect_true(any(grepl("most common log2FC is -0\\.[23]", get_run_log())))
})

test_that("TREAT at 1.5-fold is stricter than the 1.5-fold cutoff and runs end to end", {
  d <- withr_tempdir()
  cfg <- read_config(simulate_dataset(d, n_groups = 2, n_per_group = 4, seed = 11))
  r_cut <- run_pipeline(cfg, verbose = FALSE)
  cfg$statistics$significance <- "treat"
  expect_equal(cfg$statistics$treat_fc, 1.5)        # default matches the cutoff
  cfg$output$dir <- file.path(d, "treat15")
  r_t15 <- run_pipeline(cfg, verbose = FALSE)
  n <- function(r) r$de$summary$significant
  expect_lte(n(r_t15), n(r_cut))
  expect_match(significance_label(r_t15$de$settings), "TREAT, \\|FC\\| > 1.5")
  # TREAT p-values are never smaller than ordinary moderated p-values
  a <- r_cut$de$tables[[1]]; b <- r_t15$de$tables[[1]]
  m <- merge(a, b, by = "protein_id")
  ok <- !is.na(m$P.Value.x) & !is.na(m$P.Value.y)
  expect_true(all(m$P.Value.y[ok] >= m$P.Value.x[ok] - 1e-12))
  cfg$statistics$treat_fc <- 1
  expect_error(validate_config(cfg), "treat_fc")
})

test_that("detection testing finds on/off proteins and handles small groups", {
  set.seed(4)
  n <- 300
  x <- matrix(rnorm(n * 12, 22, 1), n)
  rownames(x) <- paste0("P", seq_len(n)); colnames(x) <- paste0("S", 1:12)
  x[1:10, 7:12] <- NA        # only in A (6 samples)
  x[11:15, 1:6] <- NA        # only in B
  se <- make_se(x, rep(c("A", "B"), each = 6))
  cmp <- build_comparisons(c("A", "B"), "all_vs_control", "B")
  det <- run_detection(se, cmp)
  t <- det$tables[[1]]
  expect_setequal(t$protein_id[t$significant], paste0("P", 1:15))
  expect_equal(sum(t$on_off == "only_in_A"), 10)
  expect_equal(sum(t$on_off == "only_in_B"), 5)
  expect_true(all(t$delta_fraction[t$on_off == "only_in_A"] == 1))
  # 3 vs 3: on/off is still listed but nothing can be significant
  se3 <- make_se(x[, c(1:3, 7:9)], rep(c("A", "B"), each = 3))
  det3 <- run_detection(se3, cmp)
  expect_equal(det3$summary$significant, 0L)
  expect_equal(det3$summary$min_possible_p, 0.1)
  expect_equal(det3$summary$on_off, 15L)
})

test_that("pipeline writes detection results", {
  d <- withr_tempdir()
  res <- run_pipeline(read_config(simulate_dataset(d, seed = 12)), verbose = FALSE)
  expect_false(is.null(res$detection))
  expect_true(file.exists(file.path(d, "results", "Detection_summary.csv")))
})

test_that("mode normalization removes sample offsets when nothing changes (review regression)", {
  for (cs in list(c(12, 0.5), c(5, 0.5), c(6, 0.4))) {
    set.seed(cs[1]); n <- cs[1]; off <- seq(-1, 1, length.out = n)
    x <- matrix(rnorm(3000 * n, 0, cs[2]), 3000) + rnorm(3000, 22, 2) + rep(off, each = 3000)
    rownames(x) <- paste0("P", 1:3000); colnames(x) <- paste0("S", 1:n)
    se <- make_se(x, rep(c("A", "B"), length.out = n))
    y <- SummarizedExperiment::assay(pp_normalize(se, "mode"), "norm")
    left <- matrixStats::colMedians(y - rowMeans(y))     # offsets still present
    expect_lt(max(abs(left - mean(left))), 0.12)
  }
})

test_that("detection column names are right for any group names (review regression)", {
  set.seed(6)
  x <- matrix(rnorm(60 * 6, 22), 60); rownames(x) <- paste0("P", 1:60); colnames(x) <- paste0("S", 1:6)
  x[1:5, 4:6] <- NA                                   # P1-P5 only in the first group
  for (gr in list(c("a", "b"), c("dose_a", "ctrl"), c("x_b", "y_a"))) {
    se <- make_se(x, rep(gr, each = 3))
    cmp <- build_comparisons(gr, "all_vs_control", gr[1])     # gr[2] vs gr[1]
    t <- run_detection(se, cmp)$tables[[1]]
    cols <- c(paste0(c("detected_", "n_", "fraction_", "mean_log2_"), gr[2]),
              paste0(c("detected_", "n_", "fraction_", "mean_log2_"), gr[1]))
    expect_true(all(cols %in% names(t)), info = paste(gr, collapse = "/"))
    expect_equal(anyDuplicated(names(t)), 0L)
    p1 <- t[t$protein_id == "P1", ]
    expect_equal(p1[[paste0("detected_", gr[1])]], 3)  # present in the first group
    expect_equal(p1[[paste0("detected_", gr[2])]], 0)
    expect_equal(p1$on_off, paste0("only_in_", gr[1]))
  }
})

test_that("top3 is NA where none of its three chosen precursors was seen", {
  x <- rbind(c(25, NA), c(24, NA), c(23, NA), c(NA, 15)); colnames(x) <- c("s1", "s2")
  r <- summarize_precursors(x, rep("P", 4), "top3")$x
  expect_false(is.na(r[1, "s1"])); expect_true(is.na(r[1, "s2"]))
})

test_that("centering check ignores comparisons with no tested proteins", {
  sparcproteomics:::log_reset()
  de <- list(tables = list(A_vs_B = data.frame(logFC = rep(NA_real_, 5))))
  expect_silent(sparcproteomics:::.check_centering(de, "median"))
  expect_false(any(grepl("NA in NA", get_run_log())))
})
