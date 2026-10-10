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
