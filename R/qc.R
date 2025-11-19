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
fit_sim_mean_var <- function(x, nb_overdispersion){
  x_sim <- simulate_counts(x, nb_overdispersion)
  x_sim_norm <- normalise_counts(x_sim, overdispersion = nb_overdispersion)
  x_sim_norm_cor <- correct_normalised_counts(x_sim_norm, discard_pcs)
  x_sim_norm_cor

  loess_df <- data.frame(sim.mean = rowMeans(x_sim_norm_cor), sim.var = rowVars(x_sim_norm_cor))
  loess(sim.var ~ sim.mean, loess_df, span = 0.2)
}
