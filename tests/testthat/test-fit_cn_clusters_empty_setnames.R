#tests for the empty-vector guard in fit_cn_clusters: when every
#cluster pair has zero eligible bins (very small clones, all bins
#zero-variance after bootstrap), the joint MSE array is all NaN,
#which.min() returns integer(0), and the prior code crashed at
#`setNames(numeric(0), paste0("cluster_", non_ref_clusters))` with
#"names attribute must be the same length as the vector".

#We construct a fixture where every bin has identical counts within a
#cluster (zero variance bootstrap -> is_excluded_x_y all-TRUE) and
#confirm:
#  (a) the function returns a result rather than erroring,
#  (b) the result carries scale_factors of the right length and names.

test_that("fit_cn_clusters returns NA scale_factors instead of erroring on all-empty pairs", {
  skip_if_not_installed("MatrixGenerics")

  set.seed(1)
  #3 clusters, 2 cells each: ref + 2 non-ref. Identical counts within each
  #cluster -> bootstrap means are constant -> rowVars == 0 ->
  #is_excluded_x_y is all-TRUE on every pair.
  n_bins <- 4
  clusters <- factor(c("ref", "ref", "A", "A", "B", "B"))
  is_excluded <- rep(FALSE, length(clusters))
  x <- matrix(rep(c(5, 5, 7, 7, 7, 7), each = n_bins),
              nrow = n_bins, byrow = FALSE)
  rownames(x) <- paste0("chr1.p:", seq_len(n_bins))
  colnames(x) <- paste0("c", seq_len(ncol(x)))

  res <- expect_warning(
    fit_cn_clusters(
      x, clusters, ref_cluster = "ref", is_ref_female = TRUE,
      is_excluded = is_excluded, n_bootstrap = 2,
      sf_start = 0.5, sf_end = 1.5, sf_step = 0.5,
      weight_by_precision = FALSE, weight_by_cluster = FALSE,
      fit_clusters_pairwise = TRUE
    ),
    "no finite MSE entries|no eligible bins"
  )

  expect_named(res, c("scale_factors", "mse_array"))
  expect_length(res$scale_factors, 2)
  expect_named(res$scale_factors, c("cluster_A", "cluster_B"))
  expect_true(all(is.na(res$scale_factors)))
})
