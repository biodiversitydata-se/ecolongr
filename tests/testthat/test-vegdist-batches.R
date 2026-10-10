test_that("the result does not depend on batch_rows", {
  m <- synthetic_community()
  x <- long_table(m)

  for (method in c("bray", "jaccard", "euclidean")) {
    for (binary in c(FALSE, TRUE)) {
      reference <- el_vegdist(x, method = method, binary = binary)
      for (rows in c(1, 7, 1e9)) {
        got <- el_vegdist(x, method = method, binary = binary, batch_rows = rows)
        expect_equal(as.numeric(got), as.numeric(reference), tolerance = 1e-12)
      }
    }
  }
})

test_that("batch_rows must be a positive number", {
  x <- long_table(synthetic_community())

  expect_error(el_vegdist(x, batch_rows = 0), "batch_rows")
  expect_error(el_vegdist(x, batch_rows = "a"), "batch_rows")
})

test_that("the result does not depend on threads", {
  x <- long_table(synthetic_community())

  for (method in c("bray", "jaccard", "euclidean")) {
    reference <- el_vegdist(x, method = method, threads = 1)
    expect_identical(as.numeric(el_vegdist(x, method = method, threads = 4)), as.numeric(reference))
  }
})

test_that("threads defaults to the ecolongr.threads option", {
  x <- long_table(synthetic_community())
  seen <- NULL
  original <- accumulate_pairs
  testthat::local_mocked_bindings(
    accumulate_pairs = function(..., threads) {
      seen <<- threads
      original(..., threads = threads)
    }
  )

  withr::local_options(ecolongr.threads = 2)
  el_vegdist(x)
  expect_equal(seen, 2)

  withr::local_options(ecolongr.threads = NULL)
  el_vegdist(x)
  expect_gte(seen, 1)
  expect_lte(seen, parallel::detectCores())

  el_vegdist(x, threads = 3)
  expect_equal(seen, 3)
})

test_that("threads must be a positive whole number", {
  x <- long_table(synthetic_community())

  expect_error(el_vegdist(x, threads = 0), "`threads` must be a positive whole number")
  expect_error(el_vegdist(x, threads = 1.5), "`threads` must be a positive whole number")
  expect_error(el_vegdist(x, threads = "a"), "`threads` must be a positive whole number")
})

test_that("memory_limit does not change the result", {
  x <- long_table(synthetic_community())
  reference <- el_vegdist(x)

  expect_equal(as.numeric(el_vegdist(x, memory_limit = "64MB")), as.numeric(reference))
})

test_that("memory_limit defaults to the ecolongr.memory_limit option", {
  x <- long_table(synthetic_community())
  seen <- NULL
  original <- open_taxon_stream
  testthat::local_mocked_bindings(
    open_taxon_stream = function(..., memory_limit) {
      seen <<- memory_limit
      original(..., memory_limit = memory_limit)
    }
  )

  withr::local_options(ecolongr.memory_limit = "123MB")
  el_vegdist(x)
  expect_equal(seen, "123MB")

  el_vegdist(x, memory_limit = "77MB")
  expect_equal(seen, "77MB")

  withr::local_options(ecolongr.memory_limit = NULL)
  el_vegdist(x)
  expect_match(seen, "^[0-9.]+ ?[KMG]B$")
})

test_that("memory_limit must be a single size such as '500MB'", {
  x <- long_table(synthetic_community())

  expect_error(el_vegdist(x, memory_limit = 500), "`memory_limit` must be")
  expect_error(el_vegdist(x, memory_limit = "lots"), "`memory_limit` must be")
  expect_error(el_vegdist(x, memory_limit = "1GB'; DROP TABLE v; --"), "`memory_limit` must be")
})

test_that("el_vegdist leaves the DuckDB settings of the duckplyr session alone", {
  setting <- function() {
    con <- duckplyr:::get_default_duckdb_connection()
    unlist(DBI::dbGetQuery(con, "SELECT current_setting('memory_limit') AS m, current_setting('threads') AS t"))
  }
  x <- long_table(synthetic_community())
  before <- setting()
  el_vegdist(x, memory_limit = "64MB", threads = 2)

  expect_equal(setting(), before)
})

test_that("a memory_limit far below the table size spills to disk and gives the same result", {
  set.seed(9)
  n <- 1.5e6
  x <- unique(data.frame(
    sample = sprintf("s%03d", sample.int(100, n, TRUE)),
    taxon = sprintf("t%06d", sample.int(200000, n, TRUE))
  ))
  x$quant <- rpois(nrow(x), 3) + 1
  withr::local_dir(withr::local_tempdir())

  small <- el_vegdist(x, memory_limit = "10MB", threads = 1)

  expect_equal(list.files(all.files = TRUE, no.. = TRUE), character())
  expect_equal(as.numeric(small), as.numeric(el_vegdist(x, threads = 1)))
})
