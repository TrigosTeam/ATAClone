#tests for the max_array_size guard in fit_cn_clusters. Without this
#guard, a 12-cluster sample at the default sf_step=0.05 (41 grid points)
#would allocate 41^11 ~= 5.5e17 entries (~4.4 EB) and OOM-kill R with
#no actionable error message.

test_that("fit_cn_clusters errors before allocating an oversized MSE array", {
  #Build the smallest possible inputs that pass entry-validation but
  #trigger the guard. The guard fires *before* any bootstrap work, so
  #we can use a near-empty count matrix.
  n_bins <- 5
  n_cells <- 14   #2 cells per cluster, 7 clusters
  x <- matrix(1, nrow = n_bins, ncol = n_cells)
  rownames(x) <- paste0("chr1.p:", seq_len(n_bins))
  colnames(x) <- paste0("c", seq_len(n_cells))
  clusters <- factor(rep(seq_len(7), each = 2))
  is_excluded <- rep(FALSE, n_cells)

  expect_error(
    fit_cn_clusters(
      x, clusters, ref_cluster = 1, is_ref_female = TRUE,
      is_excluded = is_excluded, n_bootstrap = 1,
      sf_start = 0.5, sf_end = 2.5, sf_step = 0.05,
      weight_by_precision = FALSE, weight_by_cluster = FALSE,
      fit_clusters_pairwise = FALSE,
      max_array_size = 1e6   #41^6 ~= 4.75e9 > 1e6 -> trips guard
    ),
    "exceeding max_array_size"
  )
})

test_that("fit_cn_clusters errors when fewer than 2 cluster levels are present", {
  x <- matrix(1, nrow = 2, ncol = 4)
  rownames(x) <- c("chr1.p:1", "chr1.p:2")
  colnames(x) <- paste0("c", 1:4)
  clusters <- factor(rep(1, 4))
  is_excluded <- rep(FALSE, 4)

  expect_error(
    fit_cn_clusters(
      x, clusters, ref_cluster = 1, is_ref_female = TRUE,
      is_excluded = is_excluded, n_bootstrap = 1,
      sf_start = 0.5, sf_end = 2.5, sf_step = 0.5,
      weight_by_precision = FALSE, weight_by_cluster = FALSE,
      fit_clusters_pairwise = FALSE
    ),
    "at least 2 cluster levels"
  )
})

test_that("fit_cn_clusters guard is silent under the budget", {
  #2 clusters at sf_step=0.5 -> 5^1 = 5 entries, well under default 1e8.
  #We don't exercise the full pipeline here (that needs realistic bin
  #counts and bootstrap variance) - just confirm the guard does not
  #pre-empt the bootstrap step with a stop().
  x <- matrix(rpois(20, 5), nrow = 5, ncol = 4)
  rownames(x) <- paste0("chr1.p:", 1:5)
  colnames(x) <- paste0("c", 1:4)
  clusters <- factor(rep(1:2, each = 2))
  is_excluded <- rep(FALSE, 4)

  #Past the guard, bootstrap / get_mse_x_ref runs. We tolerate any error
  #raised inside those (they require richer fixtures to satisfy) - the
  #point is the message is *not* the max_array_size or "at least 2" guard.
  err <- tryCatch(
    fit_cn_clusters(
      x, clusters, ref_cluster = 1, is_ref_female = TRUE,
      is_excluded = is_excluded, n_bootstrap = 1,
      sf_start = 0.5, sf_end = 2.5, sf_step = 0.5,
      weight_by_precision = FALSE, weight_by_cluster = FALSE,
      fit_clusters_pairwise = FALSE
    ),
    error = function(e) conditionMessage(e)
  )
  if (!is.null(err)) {
    expect_false(grepl("exceeding max_array_size|at least 2 cluster", err))
  }
})
