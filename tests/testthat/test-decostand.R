data(varespec, package = "vegan")

fixtures <- list(varespec = as.matrix(varespec), synthetic = synthetic_community())

for (fixture in names(fixtures)) {
  for (method in c("total", "hellinger", "log", "pa")) {
    test_that(paste0(method, " on ", fixture, " matches vegan::decostand"), {
      m <- fixtures[[fixture]]
      got <- wide_from_long(as.data.frame(el_decostand(long_table(m), method = method)), m)

      expect_equal(got, vegan_decostand(m, method), tolerance = 1e-10)
    })
  }

  test_that(paste0("rclr on ", fixture, " matches vegan::decostand on the non-zero cells"), {
    m <- fixtures[[fixture]]
    got <- as.data.frame(el_decostand(long_table(m), method = "rclr"))
    reference <- vegan_decostand(m, "rclr")

    expect_equal(wide_from_long(got, m)[m > 0], reference[m > 0], tolerance = 1e-10)
    expect_equal(nrow(got), sum(m > 0))
    expect_true(all(is.na(reference[m == 0])))
  })
}

test_that("log with another base matches vegan::decostand", {
  m <- fixtures$synthetic
  for (base in c(2, 10, exp(1))) {
    got <- wide_from_long(as.data.frame(el_decostand(long_table(m), method = "log", logbase = base)), m)
    expect_equal(got, vegan_decostand(m, "log", logbase = base), tolerance = 1e-10)
  }
})

test_that("log divides by the smallest positive value of the whole table when any value is not a whole number", {
  m <- fixtures$synthetic
  m[2, 2] <- 0
  m[5, 7] <- 0.5
  m[5, 8] <- 3.5
  got <- wide_from_long(as.data.frame(el_decostand(long_table(m), method = "log")), m)

  expect_equal(got, vegan_decostand(m, "log"), tolerance = 1e-10)
  expect_true(min(got[m > 0]) < 1)
})

test_that("log leaves whole-number data unscaled", {
  m <- fixtures$synthetic
  got <- wide_from_long(as.data.frame(el_decostand(long_table(m), method = "log")), m)

  expect_equal(got[m == 1], rep(1, sum(m == 1)))
  expect_equal(got[m == 4], rep(3, sum(m == 4)))
})

test_that("the result is a lazy duckplyr table, also for a stingy input", {
  x <- duckplyr::as_duckdb_tibble(long_table(fixtures$synthetic), prudence = "stingy")

  for (method in c("total", "hellinger", "log", "pa", "rclr")) {
    got <- el_decostand(x, method = method)
    expect_s3_class(got, "duckplyr_df")
    expect_gt(nrow(as.data.frame(got)), 0)
  }
})

test_that("the result keeps the input's column names and drops other columns", {
  x <- long_table(fixtures$synthetic)
  names(x) <- c("site", "species", "abundance")
  x$note <- "kept?"
  got <- as.data.frame(el_decostand(x, method = "total", sample = "site", taxon = "species", quant = "abundance"))

  expect_named(got, c("site", "species", "abundance"))
})

test_that("explicit zero rows are ignored and never appear in the result", {
  m <- fixtures$synthetic
  with_zeros <- rbind(long_table(m), data.frame(sample = "s04", taxon = "t59", quant = 0))

  for (method in c("total", "hellinger", "log", "pa", "rclr")) {
    got <- as.data.frame(el_decostand(with_zeros, method = method))
    expect_equal(nrow(got), sum(m > 0))
    expect_false(any(got$taxon == "t59" & got$sample == "s04"))
  }
})

test_that("total gives every sample the sum 1 and hellinger squares back to total", {
  x <- long_table(fixtures$varespec)
  total <- as.data.frame(el_decostand(x, "total"))
  hellinger <- as.data.frame(el_decostand(x, "hellinger"))

  expect_equal(unname(tapply(total$quant, total$sample, sum)), rep(1, nlevels(factor(total$sample))), tolerance = 1e-12)
  key <- function(d) paste(d$sample, d$taxon)
  expect_equal(hellinger$quant[match(key(total), key(hellinger))]^2, total$quant, tolerance = 1e-12)
})

test_that("total is idempotent", {
  once <- el_decostand(long_table(fixtures$synthetic), "total")
  twice <- el_decostand(once, "total")
  key <- function(d) paste(d$sample, d$taxon)
  a <- as.data.frame(once)
  b <- as.data.frame(twice)

  expect_equal(b$quant[match(key(a), key(b))], a$quant, tolerance = 1e-12)
})

test_that("pa is 1 on every row", {
  got <- as.data.frame(el_decostand(long_table(fixtures$synthetic), "pa"))

  expect_true(all(got$quant == 1))
})

test_that("rclr is the same for counts and relative abundances", {
  x <- long_table(fixtures$synthetic)
  relative <- as.data.frame(el_decostand(x, "total"))
  a <- as.data.frame(el_decostand(x, "rclr"))
  b <- as.data.frame(el_decostand(relative, "rclr"))
  key <- function(d) paste(d$sample, d$taxon)

  expect_equal(b$quant[match(key(a), key(b))], a$quant, tolerance = 1e-10)
})

test_that("rclr equals log(x) minus the mean log of the positive values, within each sample", {
  m <- fixtures$synthetic
  got <- wide_from_long(as.data.frame(el_decostand(long_table(m), "rclr")), m)
  independent <- t(apply(m, 1, function(r) ifelse(r > 0, log(r) - mean(log(r[r > 0])), 0)))

  expect_equal(got, independent, tolerance = 1e-10)
})

test_that("rclr is 0 for a sample with a single taxon and sums to 0 within every sample", {
  m <- fixtures$synthetic
  got <- as.data.frame(el_decostand(long_table(m), "rclr"))

  expect_equal(got$quant[got$sample == "s03"], 0)
  expect_equal(unname(tapply(got$quant, got$sample, sum)), rep(0, length(unique(got$sample))), tolerance = 1e-10)
})

test_that("decostand then vegdist gives the same distances as vegan", {
  m <- fixtures$varespec

  for (method in c("total", "hellinger")) {
    got <- el_vegdist(el_decostand(long_table(m), method), method = "euclidean")
    expect_equal(as_full_matrix(got), as_full_matrix(vegan::vegdist(vegan::decostand(m, method), "euclidean")), tolerance = 1e-10)
  }

  got <- el_vegdist(el_decostand(long_table(m), "total"), method = "bray")
  expect_equal(as_full_matrix(got), as_full_matrix(vegan::vegdist(vegan::decostand(m, "total"), "bray")), tolerance = 1e-10)
})

test_that("bad input is an error", {
  x <- long_table(fixtures$synthetic)

  expect_error(el_decostand(x, method = "chi.square"), "should be one of")
  expect_error(el_decostand(x, method = "log", logbase = 1), "`logbase` must be")
  expect_error(el_decostand(x, method = "log", logbase = "a"), "`logbase` must be")
  expect_error(el_decostand(dplyr::rename(x, abundance = quant)), "quant")
  expect_error(el_decostand(transform(x, quant = -quant)), "negative")
  expect_error(el_decostand(transform(x, quant = NA_real_)), "missing")
})

test_that("check = FALSE skips the validation query", {
  x <- duckplyr::as_duckdb_tibble(long_table(fixtures$synthetic), prudence = "stingy")
  seen <- 0
  testthat::local_mocked_bindings(
    check_quantities = function(...) {
      seen <<- seen + 1
      invisible()
    }
  )

  el_decostand(x, "total", check = FALSE)
  expect_equal(seen, 0)
  el_decostand(x, "total")
  expect_equal(seen, 1)
})
