data(varespec, package = "vegan")

fixtures <- list(
  varespec = as.matrix(varespec),
  synthetic = present_taxa(synthetic_community()),
  lmo = lmo_subset()$matrix
)

fit_matrix <- function(m, ...) {
  t <- rclr_triplets(m)
  f <- optspace_triplets(t$i, t$j, t$v, t$n, t$p, ...)
  double_centre(f$X %*% f$S %*% t(f$Y))
}

for (fixture in names(fixtures)) {
  test_that(paste0("optspace on triplets reproduces vegan's completed rclr matrix for ", fixture), {
    m <- fixtures[[fixture]]

    expect_equal(unname(fit_matrix(m)), unname(vegan_rclr_imputed(m)), tolerance = 1e-6)
  })
}

test_that("the completed matrix has rank 3 and is double-centred", {
  m <- fixtures$lmo
  completed <- fit_matrix(m)

  expect_equal(sum(svd(completed)$d > 1e-8 * svd(completed)$d[1]), 3)
  expect_equal(max(abs(rowMeans(completed))), 0, tolerance = 1e-10)
  expect_equal(max(abs(colMeans(completed))), 0, tolerance = 1e-10)
})

test_that("a different rank, iteration count and tolerance match vegan::optspace", {
  m <- fixtures$synthetic
  t <- rclr_triplets(m)
  r <- suppressWarnings(vegan::decostand(m, "rclr", impute = FALSE))

  for (cfg in list(list(r = 2, niter = 5), list(r = 4, niter = 5), list(r = 3, niter = 8))) {
    ref <- vegan::optspace(as.matrix(r), ropt = cfg$r, niter = cfg$niter, tol = 1e-5)$M
    f <- optspace_triplets(t$i, t$j, t$v, t$n, t$p, r = cfg$r, niter = cfg$niter)

    expect_equal(unname(double_centre(f$X %*% f$S %*% t(f$Y))), matrix(as.numeric(ref), nrow(ref)), tolerance = 1e-6)
  }
})

test_that("the fit does not depend on the threads or on the order of the cells", {
  m <- fixtures$lmo
  t <- rclr_triplets(m)
  reference <- optspace_triplets(t$i, t$j, t$v, t$n, t$p, threads = 1)
  shuffled <- sample(length(t$v))

  threaded <- optspace_triplets(t$i, t$j, t$v, t$n, t$p, threads = 4)
  reordered <- optspace_triplets(t$i[shuffled], t$j[shuffled], t$v[shuffled], t$n, t$p, threads = 1)
  as_m <- function(f) f$X %*% f$S %*% t(f$Y)

  expect_equal(as_m(threaded), as_m(reference), tolerance = 1e-9)
  expect_equal(as_m(reordered), as_m(reference), tolerance = 1e-9)
})

test_that("robust.aitchison matches vegan::vegdist", {
  for (m in fixtures[c("varespec", "synthetic", "lmo")]) {
    got <- el_vegdist(long_table(m), method = "robust.aitchison")
    expected <- suppressWarnings(vegan::vegdist(m, method = "robust.aitchison"))

    expect_s3_class(got, "dist")
    expect_equal(as_full_matrix(got), as_full_matrix(expected), tolerance = 1e-6)
  }
})

test_that("robust.aitchison is the Euclidean distance on the completed matrix", {
  m <- fixtures$lmo

  expect_equal(
    as_full_matrix(el_vegdist(long_table(m), method = "robust.aitchison")),
    as_full_matrix(stats::dist(vegan_rclr_imputed(m))),
    tolerance = 1e-6
  )
})

test_that("robust.aitchison does not depend on batch_rows or threads", {
  x <- long_table(fixtures$lmo)
  reference <- el_vegdist(x, method = "robust.aitchison", threads = 1)

  expect_equal(as.numeric(el_vegdist(x, method = "robust.aitchison", batch_rows = 1)), as.numeric(reference), tolerance = 1e-8)
  expect_equal(as.numeric(el_vegdist(x, method = "robust.aitchison", threads = 4)), as.numeric(reference), tolerance = 1e-8)
})

test_that("robust.aitchison rejects unsupported arguments", {
  x <- long_table(fixtures$synthetic)

  expect_error(el_vegdist(x, method = "robust.aitchison", binary = TRUE), "binary")
  expect_error(el_vegdist(x, method = "robust.aitchison", samples = c(unique(x$sample), "empty")), "all-zero samples")
  expect_error(el_vegdist(transform(x, quant = -quant), method = "robust.aitchison"), "negative")
})

test_that("el_rda with transform = 'rclr' matches vegan::rda on the completed matrix", {
  for (m in fixtures[c("varespec", "synthetic", "lmo")]) {
    got <- el_rda(long_table(m), k = 3, transform = "rclr")

    expect_equal(length(vegan::eigenvals(got)), 3)
    compare_rda(got, vegan_rclr_imputed(m), tolerance = 1e-6)
  }
})

test_that("el_rda rclr is capped at the rank of the completion and rejects scale = TRUE", {
  m <- fixtures$lmo
  got <- el_rda(long_table(m), k = 10, transform = "rclr")

  expect_equal(length(vegan::eigenvals(got)), 3)
  expect_error(el_rda(long_table(m), transform = "rclr", scale = TRUE), "scale")
  expect_error(el_rda(long_table(m), transform = "rclr", pseudocount = 1), "pseudocount")
})

test_that("rclr with imputation recovers the seasonal structure of the LMO subset as vegan's does", {
  lmo <- lmo_subset()
  meta <- utils::read.delim(testthat::test_path("fixtures", "lmo_subset_samples.tsv"), stringsAsFactors = FALSE)
  meta$doy <- as.POSIXlt(as.Date(meta$date))$yday + 1
  seasonal_r2 <- function(model) {
    sites <- vegan::scores(model, display = "sites", choices = 1:2)
    m <- meta[match(rownames(sites), meta$sample), ]
    harmonics <- cbind(sin(2 * pi * m$doy / 365), cos(2 * pi * m$doy / 365), sin(4 * pi * m$doy / 365), cos(4 * pi * m$doy / 365))
    r2 <- sapply(1:2, function(a) summary(stats::lm(sites[, a] ~ harmonics))$r.squared)
    sum(vegan::eigenvals(model)[1:2] * r2) / sum(vegan::eigenvals(model)[1:2])
  }

  got <- el_rda(lmo$long, k = 2, transform = "rclr")
  reference <- vegan::rda(vegan_rclr_imputed(lmo$matrix))

  expect_equal(seasonal_r2(got), seasonal_r2(reference), tolerance = 1e-5)
})

test_that("dimensions whose product exceeds the 32-bit integer range do not overflow", {
  set.seed(12)
  n <- 50000L
  p <- 50000L
  expect_gt(as.numeric(n) * p, .Machine$integer.max)
  i <- rep(seq_len(n) - 1L, each = 6)
  j <- sample.int(p, length(i), replace = TRUE) - 1L
  keep <- !duplicated(cbind(i, j))
  i <- i[keep]
  j <- j[keep]
  x <- matrix(rnorm(n * 2), n, 2)
  y <- matrix(rnorm(p * 2), p, 2)
  v <- rowSums(x[i + 1L, ] * y[j + 1L, ]) + rnorm(length(i), sd = 0.1)

  fit <- optspace_triplets(i, j, v, n, p)

  expect_true(all(is.finite(fit$X)))
  expect_true(all(is.finite(fit$Y)))
  expect_true(all(is.finite(fit$S)))
})

for (ropt in c(2, 5)) {
  test_that(paste0("robust.aitchison with ropt = ", ropt, " is the Euclidean distance on vegan's rank-", ropt, " completion"), {
    for (m in fixtures[c("synthetic", "lmo")]) {
      got <- el_vegdist(long_table(m), method = "robust.aitchison", ropt = ropt)

      expect_equal(as_full_matrix(got), as_full_matrix(stats::dist(vegan_rclr_imputed(m, ropt = ropt))), tolerance = 1e-6)
    }
  })

  test_that(paste0("el_rda rclr with ropt = ", ropt, " matches vegan::rda on the rank-", ropt, " completion"), {
    for (m in fixtures[c("synthetic", "lmo")]) {
      got <- el_rda(long_table(m), k = ropt, transform = "rclr", ropt = ropt)

      expect_equal(length(vegan::eigenvals(got)), ropt)
      compare_rda(got, vegan_rclr_imputed(m, ropt = ropt), k = ropt, tolerance = 1e-6)
    }
  })
}

test_that("a different ropt changes the result, and ropt = 3 is the default", {
  x <- long_table(fixtures$lmo)

  default <- el_vegdist(x, method = "robust.aitchison")
  expect_equal(as.numeric(el_vegdist(x, method = "robust.aitchison", ropt = 3)), as.numeric(default))
  expect_gt(max(abs(as.numeric(el_vegdist(x, method = "robust.aitchison", ropt = 5)) - as.numeric(default))), 1e-3)
  expect_equal(el_rda(x, k = 5, transform = "rclr", ropt = 3)$ecolongr$k, 3)
})

test_that("el_rda rclr returns at most ropt axes, however large k is", {
  got <- el_rda(long_table(fixtures$lmo), k = 50, transform = "rclr", ropt = 5)

  expect_equal(length(vegan::eigenvals(got)), 5)
})

test_that("ropt must be a whole number below the number of samples and taxa, and belongs to rclr only", {
  x <- long_table(fixtures$synthetic)

  for (bad in list(0, -1, 2.5, NA, "a", c(2, 3), 1000)) {
    expect_error(el_vegdist(x, method = "robust.aitchison", ropt = bad), "`ropt` must be")
    expect_error(el_rda(x, transform = "rclr", ropt = bad), "`ropt` must be")
  }
  expect_error(el_vegdist(x, method = "bray", ropt = 2), "ropt")
  expect_error(el_vegdist(x, method = "aitchison", ropt = 2), "ropt")
  expect_error(el_rda(x, transform = "none", ropt = 2), "ropt")
  expect_error(el_rda(x, transform = "clr", ropt = 2), "ropt")
})
