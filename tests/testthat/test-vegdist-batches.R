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
