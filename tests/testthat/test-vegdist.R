data(varespec, package = "vegan")

test_that("bray on a plain data frame matches vegan::vegdist", {
  expected <- vegan::vegdist(varespec, method = "bray")
  got <- el_vegdist(long_table(varespec), method = "bray")

  expect_s3_class(got, "dist")
  expect_equal(as_full_matrix(got), as_full_matrix(expected), tolerance = 1e-10)
})

test_that("bray on a lazy duckplyr table matches vegan::vegdist", {
  expected <- vegan::vegdist(varespec, method = "bray")
  lazy <- duckplyr::as_duckdb_tibble(long_table(varespec))
  got <- el_vegdist(lazy, method = "bray")

  expect_equal(as_full_matrix(got), as_full_matrix(expected), tolerance = 1e-10)
})

test_that("the fixture can tell counts from relative abundances", {
  on_counts <- vegan::vegdist(varespec, method = "bray")
  on_relab <- vegan::vegdist(varespec / rowSums(varespec), method = "bray")

  expect_gt(max(abs(on_counts - on_relab)), 0.01)
})

test_that("a sample listed in `samples` but absent from the table is all zero", {
  m <- rbind(varespec[1:3, ], empty = 0)
  expected <- as_full_matrix(suppressWarnings(vegan::vegdist(m, method = "bray")))
  got <- el_vegdist(long_table(varespec[1:3, ]), method = "bray", samples = rownames(m))

  expect_equal(as_full_matrix(got), expected, tolerance = 1e-10)
  expect_equal(unname(as_full_matrix(got)["empty", 1:3]), c(1, 1, 1))
})

test_that("without `samples` the all-zero sample is absent from the result", {
  got <- el_vegdist(long_table(varespec[1:3, ]), method = "bray")

  expect_equal(attr(got, "Size"), 3L)
  expect_false("empty" %in% attr(got, "Labels"))
})
