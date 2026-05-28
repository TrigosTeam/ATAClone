#tests for compute_barcode_odds: chemistry-mismatch warning + custom-table
#substitution. The bundled `barcode_probabilities` is calibrated for the
#10X Multiome v1 chemistry; users on other chemistries should either swap
#the table or accept a clear warning + NA odds.

test_that("compute_barcode_odds returns named vector with frac_known attr", {
  skip_if_not_installed("Matrix")  #lightweight dep already declared upstream
  probs <- setNames(rep(1e-6, 4), c("A", "B", "C", "D"))
  odds <- compute_barcode_odds(c("A", "B"), probabilities = probs,
                               n_total_barcodes = 1e6)
  expect_named(odds, c("A", "B"))
  expect_equal(as.numeric(odds), c(1, 1))
  expect_equal(attr(odds, "frac_known"), 1)
})

test_that("compute_barcode_odds returns NA for unknown barcodes", {
  probs <- setNames(c(2e-6, 1e-6), c("A", "B"))
  odds <- compute_barcode_odds(c("A", "Z"), probabilities = probs,
                               n_total_barcodes = 1e6)
  expect_equal(as.numeric(odds), c(2, NA))
  expect_equal(attr(odds, "frac_known"), 0.5)
})

test_that("compute_barcode_odds warns on chemistry mismatch", {
  probs <- setNames(c(1e-6, 1e-6), c("A", "B"))
  expect_warning(
    odds <- compute_barcode_odds(c("X", "Y", "Z"), probabilities = probs,
                                 n_total_barcodes = 1e6,
                                 mismatch_threshold = 0.5),
    "probability table"
  )
  expect_equal(attr(odds, "frac_known"), 0)
})

test_that("compute_barcode_odds is silent at frac_known >= threshold", {
  probs <- setNames(c(1e-6, 1e-6), c("A", "B"))
  expect_silent(
    compute_barcode_odds(c("A", "B"), probabilities = probs,
                         n_total_barcodes = 1e6,
                         mismatch_threshold = 0.5)
  )
})

test_that("compute_barcode_odds defaults to bundled multiome-v1 table", {
  data("barcode_probabilities", package = "ATAClone", envir = environment())
  some_barcodes <- names(barcode.probs)[1:5]
  odds <- compute_barcode_odds(some_barcodes)
  expect_named(odds, some_barcodes)
  expect_true(all(!is.na(odds)))
  expect_equal(attr(odds, "frac_known"), 1)
})

test_that("compute_barcode_odds errors on unnamed probabilities", {
  expect_error(
    compute_barcode_odds(c("A", "B"), probabilities = c(1e-6, 1e-6),
                         n_total_barcodes = 1e6),
    "named numeric vector"
  )
})
