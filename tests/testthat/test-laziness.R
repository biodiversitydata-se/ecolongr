data(varespec, package = "vegan")

stingy <- function(d) duckplyr::as_duckdb_tibble(d, prudence = "stingy")

test_that("duckplyr errors on a stingy frame instead of falling back to dplyr", {
  x <- stingy(data.frame(g = c("a", "a", "b"), v = 1:3))

  expect_error(dplyr::group_by(x, g), "cannot be carried out by DuckDB")
  expect_error(dplyr::mutate(x, w = v / sum(v), .by = g), "cannot be carried out by DuckDB")
  expect_no_error(dplyr::mutate(x, w = v / sum(v, na.rm = TRUE), .by = g))
})

test_that("el_vegdist stays inside DuckDB on a stingy frame", {
  expected <- vegan::vegdist(varespec, method = "bray")
  got <- el_vegdist(stingy(long_table(varespec)), method = "bray")

  expect_equal(as_full_matrix(got), as_full_matrix(expected), tolerance = 1e-10)
})

test_that("el_vegdist stays inside DuckDB when `samples` is given", {
  got <- el_vegdist(
    stingy(long_table(varespec[1:3, ])),
    samples = c(rownames(varespec)[1:3], "empty")
  )

  expect_equal(attr(got, "Size"), 4L)
})

test_that("el_vegdist converts its input with prudence 'stingy'", {
  seen <- NULL
  original <- duckplyr::as_duckdb_tibble
  testthat::local_mocked_bindings(
    as_duckdb_tibble = function(x, ..., prudence = "lavish") {
      seen <<- prudence
      original(x, ..., prudence = prudence)
    },
    .package = "duckplyr"
  )

  el_vegdist(long_table(varespec[1:3, ]))

  expect_equal(seen, "stingy")
})
