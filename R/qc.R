#' @export
get_expected_non_zeros <- function(x, overdispersion){
  lib.sizes <- colSums(x)
  gene.props <- rowSums(x)/sum(x)
  mu.mat <- t(lib.sizes %*% t(gene.props))
  mu.list <- as.list(as.data.frame(mu.mat))
  expected.zeros <- numeric()
  if (overdispersion == 0){
    for (i in seq_along(mu.list)){
      expected.zeros[i] <- sum(sapply(mu.list[[i]], ppois, q = 0))
    }
  }
  else {
    for (i in seq_along(mu.list)){
      expected.zeros[[i]] <- sum(sapply((1 / overdispersion)/(mu.list[[i]] +  1 / overdispersion),
                                        pnbinom, q = 0, size = 1 / overdispersion))
    }
  }
  nrow(x) - expected.zeros
}

#this is for fitting the same bin across cells
get_glm <- function(y, x, y.offset){
  y.glm <- glm(y ~ log(x) + offset(log(y.offset)), family = "poisson")
  y.glm
}

get_pred <- function(y.glm, x){
  y.pred <- exp(y.glm$coefficients[1] + y.glm$coefficients[2] * log(x))
  y.pred
}

#' @export
fit_stable_frac_regression <- function(y, stable_counts, all_counts){
  x = colSums(stable_counts) / colSums(all_counts)
  y.list <- as.list(as.data.frame(t(as.matrix(y))))
  y.offset <- colSums(y)
  y.glm.list <- lapply(y.list, get_glm, x = x, y.offset = y.offset)
  y.pred.list <- lapply(y.glm.list, get_pred, x = x)
  y.pred.df <- reshape2::melt(y.pred.list)
  colnames(y.pred.df) <- c("pred_bin_frac","bin")
  y.pred.df$cell_stable_frac <- x
  y.pred.df$bin_stable_frac <- rep(rowSums(stable_counts) / rowSums(all_counts), each = ncol(y))
  y.pred.df
}

#' @export
fit_sim_mean_var <- function(x, nb_overdispersion, pseudo_count, lambda){
  x_sim <- simulate_counts(x, nb_overdispersion)
  x_sim_norm <- normalise_counts(x_sim, overdispersion = nb_overdispersion, pseudo_count, lambda)
  x_sim_norm_cor <- correct_normalised_counts(x_sim_norm, 1)
  x_sim_norm_cor

  loess_df <- data.frame(sim.mean = rowMeans(x_sim_norm_cor), sim.var = rowVars(x_sim_norm_cor))
  loess(sim.var ~ sim.mean, loess_df, span = 0.2)
}

#this is for fitting across different bins in the same cell. Example inputs:
#x <- log(rowMeans(stable_counts_ref)[rowMeans(stable_counts_ref) > 0])
#y <- as.list(as.data.frame(as.matrix(stable_counts_filtered)[rowMeans(stable_counts_ref) > 0,]))
#or
#bin_fracs <- rowSums(stable_counts_filtered) / rowSums(all_counts_filtered)
#x <- log(bin_fracs[bin_fracs > 0])
#y <- as.list(as.data.frame(as.matrix(stable_counts_filtered)[bin_fracs > 0,]))
get_pois_reg <- function(y, x){
  glm(y ~ x, poisson(link = "log"))$coefficients
}

#' @export
get_pois_regression_betas <- function(y, x){
  do.call(rbind, lapply(y, get_pois_reg, x))[,2]
}

#' Compute per-barcode odds against a chemistry-specific probability table
#'
#' Returns the ratio of each barcode's empirical probability (under the
#' supplied table) to the uniform prior over the chemistry's whitelist size.
#' The bundled `barcode_probabilities` table is calibrated against the 10X
#' Multiome v1 chemistry; users on a different chemistry (standalone
#' CellRanger-ATAC, scATAC-seq via a different platform) should pass their
#' own `probabilities` table or expect a warning + NA odds.
#'
#' @param barcodes character vector of cell barcodes to score.
#' @param probabilities a named numeric vector mapping barcode -> probability.
#'   When NULL (the default), the bundled multiome v1 table is loaded.
#' @param n_total_barcodes the chemistry's whitelist size (the denominator of
#'   the uniform prior). Defaults to 735320 (10X Multiome v1 ATAC whitelist).
#' @param mismatch_threshold if the fraction of `barcodes` covered by
#'   `probabilities` falls below this value, a warning is emitted that the
#'   table likely doesn't match the assay's chemistry. Default 0.5.
#' @return a named numeric vector of barcode odds, with NA for any barcode
#'   not present in the probability table. `attr(., "frac_known")` is set
#'   so callers can decide whether to apply the `barcode_odds` filter.
#' @export
compute_barcode_odds <- function(barcodes, probabilities = NULL,
                                 n_total_barcodes = 735320,
                                 mismatch_threshold = 0.5) {
  if (is.null(probabilities)) {
    e <- new.env()
    utils::data("barcode_probabilities", package = "ATAClone", envir = e)
    probabilities <- e$barcode.probs
  }
  if (is.null(names(probabilities))) {
    stop("`probabilities` must be a named numeric vector (barcode -> prob).")
  }
  odds <- n_total_barcodes * probabilities[barcodes]
  names(odds) <- barcodes
  frac_known <- mean(!is.na(odds))
  attr(odds, "frac_known") <- frac_known
  if (frac_known < mismatch_threshold) {
    warning(sprintf(
      paste0("Only %.1f%% of barcodes found in the probability table ",
             "(threshold %.1f%%). The bundled `barcode_probabilities` is ",
             "calibrated for 10X Multiome v1; if your data uses a different ",
             "chemistry, supply a matching `probabilities` table or skip ",
             "the barcode-odds filter."),
      100 * frac_known, 100 * mismatch_threshold
    ))
  }
  odds
}
