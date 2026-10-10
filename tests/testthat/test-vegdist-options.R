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

signed_community <- function() {
  set.seed(8)
  m <- synthetic_community()
  nz <- which(m > 0)
  m[nz] <- m[nz] * sample(c(-1, 1), length(nz), TRUE) + sample(c(0.25, -0.5), length(nz), TRUE)
  m
}

test_that("euclidean accepts negative quantities and matches vegan", {
  m <- signed_community()
  expect_true(any(m < 0))

  got <- el_vegdist(long_table_nonzero(m), method = "euclidean")

  expect_equal(as_full_matrix(got), as_full_matrix(vegan::vegdist(m, "euclidean")), tolerance = 1e-10)
})

test_that("negative quantities are still an error for bray, jaccard and binary euclidean", {
  x <- long_table_nonzero(signed_community())

  expect_error(el_vegdist(x, method = "bray"), "negative")
  expect_error(el_vegdist(x, method = "jaccard"), "negative")
  expect_error(el_vegdist(x, method = "euclidean", binary = TRUE), "negative")
})

test_that("negative quantities and zero rows together: zeros are still ignored for euclidean", {
  m <- signed_community()
  x <- rbind(long_table_nonzero(m), data.frame(sample = "s04", taxon = "t59", quant = 0))

  expect_equal(
    as_full_matrix(el_vegdist(x, method = "euclidean")),
    as_full_matrix(vegan::vegdist(m, "euclidean")),
    tolerance = 1e-10
  )
})

test_that("rclr then euclidean is vegan's robust Aitchison distance when there are no zeros", {
  set.seed(4)
  m <- matrix(rpois(15 * 12, 6) + 1, 15, 12, dimnames = list(sprintf("s%02d", 1:15), sprintf("t%02d", 1:12)))

  got <- el_vegdist(el_decostand(long_table(m), "rclr"), method = "euclidean")

  expect_equal(as_full_matrix(got), as_full_matrix(vegan::vegdist(m, "robust.aitchison")), tolerance = 1e-10)
})

test_that("rclr then euclidean treats absent cells as 0 when there are zeros", {
  for (m in list(as.matrix(varespec), synthetic_community())) {
    rclr <- vegan_decostand(m, "rclr")
    rclr[is.na(rclr)] <- 0

    got <- el_vegdist(el_decostand(long_table(m), "rclr"), method = "euclidean")

    expect_equal(as_full_matrix(got), as_full_matrix(vegan::vegdist(rclr, "euclidean")), tolerance = 1e-10)
  }
})

test_that("a sample that appears only with zero rows is kept as an all-zero sample", {
  m <- synthetic_community()[1:4, ]
  x <- rbind(long_table(m), data.frame(sample = "zero", taxon = c("t01", "t02"), quant = 0))
  with_row <- rbind(m, zero = 0)

  for (method in c("bray", "jaccard", "euclidean")) {
    got <- el_vegdist(x, method = method)
    expect_equal(
      as_full_matrix(got),
      as_full_matrix(suppressWarnings(vegan::vegdist(with_row, method))),
      tolerance = 1e-10
    )
  }
})

test_that("integer counts above the 32-bit square limit do not overflow", {
  set.seed(2)
  m <- matrix(0L, 12, 20, dimnames = list(sprintf("s%02d", 1:12), sprintf("t%02d", 1:20)))
  m[sample(length(m), 90)] <- sample(c(60000L, 120000L, 3L, 90000L), 90, TRUE)
  x <- long_table_nonzero(m)
  x$quant <- as.integer(x$quant)
  expect_type(x$quant, "integer")

  for (method in c("euclidean", "aitchison", "bray", "jaccard")) {
    expected <- if (method == "aitchison") vegan::vegdist(m[, colSums(m) > 0], "aitchison", pseudocount = 1) else vegan::vegdist(m, method)
    expect_equal(as_full_matrix(el_vegdist(x, method = method)), as_full_matrix(expected), tolerance = 1e-10)
  }

  rda <- el_rda(x, k = 3, solver = "exact")
  compare_rda(rda, m, k = 3)
  rda_clr <- el_rda(x, k = 3, transform = "clr", solver = "exact")
  compare_rda(rda_clr, vegan::decostand(m[, colSums(m) > 0], "clr", pseudocount = 1))
})
