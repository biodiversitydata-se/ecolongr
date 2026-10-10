batch_table <- function() {
  data.frame(
    sample = c("a", "b", "c", "a", "b", "a", "c", "d", "e", "e"),
    taxon  = c("t3", "t3", "t3", "t1", "t1", "t2", "t2", "t2", "t2", "t4"),
    quant  = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  )
}

collect_batches <- function(x, batch_rows) {
  out <- list()
  for_each_taxon_batch(x, batch_rows, function(b) out[[length(out) + 1]] <<- b)
  out
}

test_that("the batches together are the table sorted by taxon", {
  x <- batch_table()
  got <- do.call(rbind, collect_batches(x, batch_rows = 3))

  expect_equal(sort(paste(got$taxon, got$sample)), sort(paste(x$taxon, x$sample)))
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
