#tests for read_matrix's tolerance of both standard 10X MTX (no header
#row in barcodes/features.tsv.gz) and the legacy in-package sentinel
#header ("colnames.x." / "rownames.x.").

write_mtx_triplet <- function(dir, mat, barcodes, features) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  Matrix::writeMM(methods::as(mat, "dgTMatrix"), file.path(dir, "matrix.mtx"))
  R.utils::gzip(file.path(dir, "matrix.mtx"), overwrite = TRUE)
  writeLines(barcodes, gzfile(file.path(dir, "barcodes.tsv.gz")))
  writeLines(features, gzfile(file.path(dir, "features.tsv.gz")))
}

test_that("read_matrix handles a standard 10X MTX (no sentinel header)", {
  skip_if_not_installed("R.utils")

  tmp <- tempfile("mtx_no_header_")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  barcodes <- c("AAAA-1", "CCCC-1", "GGGG-1")
  features <- c("chr1.p:1-1e6", "chr1.p:1e6-2e6")
  m <- Matrix::sparseMatrix(i = c(1, 1, 1), j = c(1, 2, 3), x = c(1, 2, 3), dims = c(2, 3))

  write_mtx_triplet(tmp, m, barcodes, features)

  x <- read_matrix(tmp)
  expect_equal(colnames(x), barcodes)
  expect_equal(rownames(x), features)
  expect_equal(dim(x), c(length(features), length(barcodes)))
})

test_that("read_matrix tolerates the legacy sentinel header", {
  skip_if_not_installed("R.utils")

  tmp <- tempfile("mtx_sentinel_")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  barcodes <- c("AAAA-1", "CCCC-1")
  features <- c("chr1.p:1-1e6", "chr1.p:1e6-2e6")
  m <- Matrix::sparseMatrix(i = c(1, 2, 1), j = c(1, 2, 2), x = c(1, 2, 3), dims = c(2, 2))

  write_mtx_triplet(tmp, m,
                   barcodes = c("colnames.x.", barcodes),
                   features = c("rownames.x.", features))

  x <- read_matrix(tmp)
  expect_equal(colnames(x), barcodes)
  expect_equal(rownames(x), features)
})

test_that("read_matrix errors on a name/dim mismatch", {
  skip_if_not_installed("R.utils")

  tmp <- tempfile("mtx_mismatch_")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  barcodes <- c("AAAA-1", "CCCC-1", "GGGG-1", "TTTT-1")  #1 too many
  features <- c("chr1.p:1-1e6", "chr1.p:1e6-2e6")
  m <- Matrix::sparseMatrix(i = c(1, 1, 1), j = c(1, 2, 3), x = c(1, 2, 3), dims = c(2, 3))

  write_mtx_triplet(tmp, m, barcodes, features)

  expect_error(read_matrix(tmp), "does not match matrix dimension")
})
