lv <- c("Control", "Low", "High")

test_that("all_vs_control compares every group with the control", {
  cmp <- build_comparisons(lv, "all_vs_control")
  expect_equal(cmp$name, c("Low_vs_Control", "High_vs_Control"))
  expect_true(all(cmp$denominator == "Control"))
})

test_that("pairwise gives k(k-1)/2 comparisons", {
  expect_equal(nrow(build_comparisons(c(lv, "Max"), "pairwise")), 6)
})

test_that("custom comparisons accept 'A vs B' and 'A - B'", {
  cmp <- build_comparisons(lv, list("High vs Low", "Low - Control"))
  expect_equal(cmp$numerator, c("High", "Low"))
  expect_equal(cmp$denominator, c("Low", "Control"))
  expect_error(build_comparisons(lv, "High vs Medium"), "unknown group")
  expect_error(build_comparisons(lv, "High versus Low"), "Cannot parse")
})

test_that("sample matching tolerates paths, extensions and Spectronaut prefixes", {
  cols <- c("D:/runs/S01.raw", "[2] S02.d", "s03.WIFF", "Unknown")
  expect_equal(match_samples(cols, c("S03", "S02", "S01")), c(3L, 2L, 1L, NA))
})
