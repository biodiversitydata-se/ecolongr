data(varespec, package = "vegan")

fixtures <- list(varespec = as.matrix(varespec), synthetic = synthetic_community())
methods <- c("bray", "jaccard", "euclidean")

for (fixture in names(fixtures)) {
  for (method in methods) {
    for (binary in c(FALSE, TRUE)) {
      test_that(paste0(method, ", binary = ", binary, " on ", fixture, " matches vegan::vegdist"), {
        m <- fixtures[[fixture]]
        expected <- vegan_reference(m, method, binary)
        got <- el_vegdist(long_table(m), method = method, binary = binary)

        expect_equal(as_full_matrix(got), as_full_matrix(expected), tolerance = 1e-10)
      })
    }
  }
}

for (method in methods) {
  for (binary in c(FALSE, TRUE)) {
    test_that(paste0(method, ", binary = ", binary, ": two all-zero samples listed in `samples` match vegan"), {
      m <- synthetic_community()
      empty <- c("empty1", "empty2")
      expected <- vegan_reference(m, method, binary, empty)
      got <- el_vegdist(long_table(m), method = method, binary = binary, samples = c(rownames(m), empty))

      expect_equal(as_full_matrix(got), as_full_matrix(expected), tolerance = 1e-10)
    })
  }
}

test_that("two all-zero samples are NaN for bray and jaccard but 0 for euclidean", {
  m <- synthetic_community()[1:3, ]
  samples <- c(rownames(m), "empty1", "empty2")

  for (method in c("bray", "jaccard")) {
    got <- as.matrix(el_vegdist(long_table(m), method = method, samples = samples))
    expect_true(is.nan(got["empty1", "empty2"]))
    expect_equal(unname(got["empty1", rownames(m)]), c(1, 1, 1))
  }

  got <- as.matrix(el_vegdist(long_table(m), method = "euclidean", samples = samples))
  expect_equal(got["empty1", "empty2"], 0)
})

test_that("euclidean between an all-zero sample and another is the norm of the other", {
  m <- synthetic_community()[1:3, ]
  got <- as.matrix(el_vegdist(long_table(m), method = "euclidean", samples = c(rownames(m), "empty")))

  expect_equal(unname(got["empty", rownames(m)]), unname(sqrt(rowSums(m^2))), tolerance = 1e-10)
})

test_that("identical samples are at distance 0", {
  m <- synthetic_community()
  for (method in methods) {
    got <- as.matrix(el_vegdist(long_table(m), method = method))
    expect_equal(got["s01", "s02"], 0)
  }
})

test_that("samples sharing no taxon are at distance 1 for bray and jaccard", {
  m <- rbind(a = c(1, 2, 0, 0), b = c(0, 0, 3, 4))
  colnames(m) <- paste0("t", 1:4)

  for (method in c("bray", "jaccard")) {
    expect_equal(as.numeric(el_vegdist(long_table(m), method = method)), 1)
  }
})

test_that("explicit zero rows in the table are ignored", {
  m <- synthetic_community()
  with_zeros <- rbind(long_table(m), data.frame(sample = "s04", taxon = "t59", quant = 0))

  for (method in methods) {
    for (binary in c(FALSE, TRUE)) {
      expect_equal(
        as_full_matrix(el_vegdist(with_zeros, method = method, binary = binary)),
        as_full_matrix(el_vegdist(long_table(m), method = method, binary = binary))
      )
    }
  }
})

test_that("the sample order follows `samples` when given", {
  m <- synthetic_community()[1:4, ]
  got <- el_vegdist(long_table(m), samples = rev(rownames(m)))

  expect_equal(attr(got, "Labels"), rev(rownames(m)))
})

test_that("the sample order is alphabetical when `samples` is not given", {
  m <- synthetic_community()[4:1, ]
  got <- el_vegdist(long_table(m))

  expect_equal(attr(got, "Labels"), sort(rownames(m)))
})

test_that("bad input is an error", {
  x <- long_table(synthetic_community())

  expect_error(el_vegdist(x, method = "manhattan"), "should be one of")
  expect_error(el_vegdist(x, samples = c("s01", "s02")), "not in `samples`")
  expect_error(el_vegdist(dplyr::rename(x, abundance = quant)), "quant")
  expect_error(el_vegdist(transform(x, quant = -quant)), "negative")
  expect_error(el_vegdist(transform(x, quant = NA_real_)), "missing")
})
