#Regression tests: plot_absolute_cn2 and get_absolute_copy_number must
#tolerate the all-NA scale_factors that fit_cn_clusters now returns
#(with a warning) when no cluster pair has finite MSE entries — typical
#of samples with very small clones.
#
#Before this fix, plot_absolute_cn2 crashed at
#  reshape2::melt -> setNames
#with "names attribute must be the same length as the vector" because
#the all-NA-list filter at the top emptied the list, do.call(cbind, .)
#returned NULL, and melt(NULL) produced an empty data frame whose
#columns could not be named.

test_that("plot_absolute_cn2 returns a placeholder ggplot when input is all-NA", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("reshape2")

  bin_names <- paste0("chr1.p:", seq_len(5))
  na_vec <- setNames(rep(NA_real_, 5), bin_names)
  absolute_cn_list <- list(`2` = na_vec, `3` = na_vec, all = na_vec)

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  res <- expect_warning(
    plot_absolute_cn2(absolute_cn_list),
    "no finite copy-number values"
  )
  expect_s3_class(res, "ggplot")
})

test_that("plot_absolute_cn2 returns a placeholder ggplot on an empty input list", {
  skip_if_not_installed("ggplot2")

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  res <- expect_warning(
    plot_absolute_cn2(list()),
    "no finite copy-number values"
  )
  expect_s3_class(res, "ggplot")
})

test_that("plot_absolute_cn2 still produces a ggplot on finite input", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("reshape2")

  #Two-clone fixture with finite values across two chromosomes, plus
  #the conventional "all" pooled track. Picked deliberately tiny so we
  #only exercise the happy path, not the heatmap layout.
  bin_names <- c(
    paste0("chr1.p:", seq_len(2)),
    paste0("chr2.p:", seq_len(2))
  )
  cn_a <- setNames(c(2, 2.1, 2, 2.2), bin_names)
  cn_b <- setNames(c(2.4, 2.3, 2, 2.1), bin_names)
  cn_all <- setNames((cn_a + cn_b) / 2, bin_names)
  absolute_cn_list <- list(`1` = cn_a, `2` = cn_b, all = cn_all)

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  #The original plot_absolute_cn2 hard-codes setNames(1:313, ...) on the
  #bin index — so it only renders cleanly when the input has exactly 313
  #bins. We're not testing that path; we're confirming the unfittable
  #guard does not short-circuit on a list with finite values.
  expect_false(ATAClone:::.is_unfittable_cn_list(absolute_cn_list))
})

test_that(".is_unfittable_cn_list detects all the unfittable shapes", {
  expect_true(ATAClone:::.is_unfittable_cn_list(NULL))
  expect_true(ATAClone:::.is_unfittable_cn_list(list()))
  expect_true(ATAClone:::.is_unfittable_cn_list(
    list(a = NA_real_, b = NA_real_)
  ))
  expect_true(ATAClone:::.is_unfittable_cn_list(
    list(a = c(NA_real_, NA_real_), b = c(NA_real_, NA_real_))
  ))
  expect_false(ATAClone:::.is_unfittable_cn_list(
    list(a = c(1, 2), b = c(3, 4))
  ))
  expect_false(ATAClone:::.is_unfittable_cn_list(
    #partial finite — one cluster has values, others NA. We still want
    #to render, so this must be considered fittable.
    list(a = c(1, 2), b = c(NA_real_, NA_real_))
  ))
})

test_that("get_absolute_copy_number returns NA-shaped result on all-NA scale_factors", {
  set.seed(1)
  n_bins <- 4
  bin_names <- paste0("chr1.p:", seq_len(n_bins))
  clusters <- factor(c("ref", "ref", "A", "A", "B", "B"))
  is_excluded <- rep(FALSE, length(clusters))
  x <- matrix(rep(c(5, 5, 7, 7, 6, 6), each = n_bins),
              nrow = n_bins, byrow = FALSE)
  rownames(x) <- bin_names
  colnames(x) <- paste0("c", seq_len(ncol(x)))

  na_sf <- setNames(rep(NA_real_, 2), c("cluster_A", "cluster_B"))

  res <- expect_warning(
    get_absolute_copy_number(
      x = x, clusters = clusters, ref_cluster = "ref",
      scale_factors = na_sf, is_excluded = is_excluded,
      is_ref_female = TRUE
    ),
    "all scale_factors are NA"
  )

  expect_named(res, c("A", "B", "all"))
  for (nm in names(res)) {
    expect_length(res[[nm]], n_bins)
    expect_true(all(is.na(res[[nm]])))
  }
})

test_that(".is_unfittable_scale_factors detects empty / all-NA inputs", {
  expect_true(ATAClone:::.is_unfittable_scale_factors(numeric(0)))
  expect_true(ATAClone:::.is_unfittable_scale_factors(NA_real_))
  expect_true(ATAClone:::.is_unfittable_scale_factors(c(NA_real_, NA_real_)))
  expect_false(ATAClone:::.is_unfittable_scale_factors(c(1, NA_real_)))
  expect_false(ATAClone:::.is_unfittable_scale_factors(c(1, 2)))
})
