collect_batches <- function(x, batch_rows, samples = NULL) {
  out <- list()
  for_each_taxon_batch(x, batch_rows, function(b) out[[length(out) + 1]] <<- b, samples = samples)
  out
}

test_that("a batch has integer sample ids, integer taxon ids and a quantity", {
  b <- collect_batches(batch_table(), batch_rows = 100)[[1]]

  expect_named(b, c("sample", "taxon", "quant"))
  expect_type(b$sample, "integer")
  expect_type(b$taxon, "integer")
  expect_type(b$quant, "double")
})

test_that("sample ids index `samples`, which defaults to the sorted sample names", {
  x <- batch_table()
  got <- do.call(rbind, collect_batches(x, batch_rows = 100))
  default <- sort(unique(x$sample))

  expect_equal(sort(default[got$sample + 1L]), sort(x$sample))
  expect_equal(min(got$sample), 0L)
  expect_equal(max(got$sample), length(default) - 1L)

  reversed <- do.call(rbind, collect_batches(x, batch_rows = 100, samples = rev(default)))
  expect_equal(sort(rev(default)[reversed$sample + 1L]), sort(x$sample))
  expect_false(identical(sort(reversed$sample), sort(got$sample)) && identical(reversed$sample, got$sample))
})

test_that("the batches together hold every row, with one taxon id per taxon", {
  x <- batch_table()
  got <- do.call(rbind, collect_batches(x, batch_rows = 3))

  expect_equal(nrow(got), nrow(x))
  expect_equal(sort(unname(table(got$taxon))), sort(unname(table(x$taxon))))
  expect_equal(sum(got$quant), sum(x$quant))
  expect_false(is.unsorted(got$taxon))
})

test_that("no taxon is split over two batches", {
  batches <- collect_batches(batch_table(), batch_rows = 3)
  taxa <- lapply(batches, function(b) unique(b$taxon))

  expect_equal(anyDuplicated(unlist(taxa)), 0L)
})

test_that("a batch holds at most batch_rows rows unless one taxon is larger", {
  batches <- collect_batches(batch_table(), batch_rows = 3)
  sizes <- vapply(batches, nrow, integer(1))
  largest_taxon <- max(table(batch_table()$taxon))

  expect_gt(length(batches), 1)
  expect_true(all(sizes <= max(3, largest_taxon)))
})

test_that("batch_rows = 1 gives one batch per taxon", {
  batches <- collect_batches(batch_table(), batch_rows = 1)

  expect_equal(length(batches), 4)
  expect_equal(sort(vapply(batches, nrow, integer(1))), c(1, 2, 3, 4))
})

test_that("a large batch_rows gives a single batch", {
  expect_equal(length(collect_batches(batch_table(), batch_rows = 1e6)), 1)
})

test_that("an empty table gives no batches", {
  expect_equal(length(collect_batches(batch_table()[0, ], batch_rows = 3)), 0)
})

test_that("a lazy stingy frame is accepted", {
  x <- duckplyr::as_duckdb_tibble(batch_table(), prudence = "stingy")

  expect_equal(sum(vapply(collect_batches(x, batch_rows = 3), nrow, integer(1))), 10)
})

test_that("a stream leaves no Parquet files behind", {
  before <- list.files(tempdir(), pattern = "\\.parquet$", recursive = TRUE)
  collect_batches(batch_table(), batch_rows = 3)
  el_vegdist(long_table(synthetic_community()))

  expect_equal(list.files(tempdir(), pattern = "\\.parquet$", recursive = TRUE), before)
})
