#' Per-barcode probabilities for the 10X Multiome v1 chemistry
#'
#' A named numeric vector of empirically-derived per-barcode probabilities
#' for the 10X Multiome v1 ATAC barcode whitelist (`barcode.probs`, length
#' 736320). Used by [compute_barcode_odds()] to derive the
#' `barcode_odds` cell-quality metric: the ratio of a barcode's probability
#' to the uniform prior `1 / 735320`.
#'
#' These probabilities are **chemistry-specific**: they encode positional
#' biases of the Multiome v1 i7/i5 indexing scheme. Users running ATAClone
#' on data from a different chemistry — standalone CellRanger-ATAC, or any
#' non-10X-Multiome-v1 platform — will see most of their barcodes outside
#' this table, producing NA odds. In that case either:
#'
#' - skip the `barcode_odds` filter, or
#' - supply a chemistry-matched probability table to
#'   [compute_barcode_odds()] via the `probabilities` argument.
#'
#' [compute_barcode_odds()] emits a warning when the fraction of supplied
#' barcodes found in the table falls below a configurable threshold.
#'
#' @docType data
#' @keywords datasets
#' @name barcode_probabilities
#' @usage data("barcode_probabilities")
#' @format A named numeric vector of length 736320; names are 16-bp barcode
#'   sequences with the standard `-1` suffix; values are probabilities
#'   summing to ~1.
#' @source Empirically estimated from a multiome v1 dataset; see the
#'   package vignette for the calibration procedure.
NULL
