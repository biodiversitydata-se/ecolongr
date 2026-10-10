data(varespec, package = "vegan")

fixtures <- list(
  varespec = as.matrix(varespec),
  synthetic = present_taxa(synthetic_community()),
  lmo = lmo_subset()$matrix
)

for (fixture in names(fixtures)) {
  for (pc in c(1, 0.5, 0.01)) {
    test_that(paste0("aitchison with pseudocount ", pc, " on ", fixture, " matches vegan::vegdist"), {
      m <- fixtures[[fixture]]
      got <- el_vegdist(long_table(m), method = "aitchison", pseudocount = pc)
      expected <- vegan::vegdist(m, method = "aitchison", pseudocount = pc)

      expect_s3_class(got, "dist")
      expect_equal(as_full_matrix(got), as_full_matrix(expected), tolerance = 1e-10)
    })
  }
}

test_that("aitchison is Euclidean distance on the clr matrix, computed independently", {
  m <- fixtures$lmo
  got <- el_vegdist(long_table(m), method = "aitchison", pseudocount = 0.1)

  expect_equal(as_full_matrix(got), as_full_matrix(stats::dist(clr_reference(m, 0.1))), tolerance = 1e-10)
})

test_that("the default pseudocount is 1", {
  m <- fixtures$synthetic

  expect_equal(
    as.numeric(el_vegdist(long_table(m), method = "aitchison")),
    as.numeric(el_vegdist(long_table(m), method = "aitchison", pseudocount = 1))
  )
})

test_that("a named pseudocount per sample matches the independent clr", {
  m <- fixtures$lmo
  pc <- stats::setNames(1 / rowSums(m), rownames(m))
  got <- el_vegdist(long_table(m), method = "aitchison", pseudocount = pc)

  expect_equal(as_full_matrix(got), as_full_matrix(stats::dist(clr_reference(m, pc))), tolerance = 1e-10)
})

test_that("pseudocount = 'inverse_depth' is 1/N_i per sample, in count space", {
  m <- fixtures$lmo
  pc <- stats::setNames(1 / rowSums(m), rownames(m))

  by_name <- el_vegdist(long_table(m), method = "aitchison", pseudocount = pc)
  shorthand <- el_vegdist(long_table(m), method = "aitchison", pseudocount = "inverse_depth")

  expect_equal(as.numeric(shorthand), as.numeric(by_name), tolerance = 1e-12)
  expect_equal(as_full_matrix(shorthand), as_full_matrix(stats::dist(clr_reference(m, pc))), tolerance = 1e-10)
})

test_that("inverse_depth makes a shallow sample's pseudocount larger than a deep one's", {
  m <- fixtures$lmo
  shallow <- which.min(rowSums(m))
  deep <- which.max(rowSums(m))
  pc <- 1 / rowSums(m)

  expect_gt(pc[shallow], pc[deep])
  expect_equal(unname(rowSums(m)[shallow] * pc[shallow]), 1)
})

test_that("an all-zero sample under inverse_depth does not produce an infinite pseudocount", {
  m <- present_taxa(fixtures$synthetic[1:5, ])
  got <- el_vegdist(long_table(m), method = "aitchison", pseudocount = "inverse_depth", samples = c(rownames(m), "empty"))

  expect_true(all(is.finite(as.matrix(got))))
})

test_that("a pseudocount that depends on the sample changes the distances", {
  m <- fixtures$lmo
  per_sample <- el_vegdist(long_table(m), method = "aitchison", pseudocount = stats::setNames(1 / rowSums(m), rownames(m)))
  scalar <- el_vegdist(long_table(m), method = "aitchison", pseudocount = 1)

  expect_gt(max(abs(per_sample - scalar)), 0.1)
})

test_that("all-zero samples listed in `samples` sit at the origin and match vegan", {
  m <- fixtures$synthetic[1:6, ]
  m <- m[, colSums(m) > 0]
  with_zero <- rbind(m, empty1 = 0, empty2 = 0)
  got <- el_vegdist(long_table(m), method = "aitchison", pseudocount = 1, samples = rownames(with_zero))
  expected <- vegan::vegdist(with_zero, method = "aitchison", pseudocount = 1)

  expect_equal(as_full_matrix(got), as_full_matrix(expected), tolerance = 1e-10)
  expect_equal(as.matrix(got)["empty1", "empty2"], 0)
})

test_that("the taxa of the table are the taxa of the clr: a taxon in no sample is not counted", {
  m <- fixtures$synthetic
  with_empty_taxon <- cbind(m, t_unseen = 0)
  without <- vegan::vegdist(m, "aitchison", pseudocount = 1)
  with <- vegan::vegdist(with_empty_taxon, "aitchison", pseudocount = 1)
  got <- el_vegdist(long_table(m), method = "aitchison", pseudocount = 1)

  expect_equal(as_full_matrix(got), as_full_matrix(without), tolerance = 1e-10)
  expect_gt(max(abs(without - with)), 1e-6)
})

test_that("the result does not depend on batch_rows or threads", {
  x <- long_table(fixtures$lmo)
  reference <- el_vegdist(x, method = "aitchison", pseudocount = 0.5, threads = 1)

  for (rows in c(1, 500)) {
    expect_equal(as.numeric(el_vegdist(x, method = "aitchison", pseudocount = 0.5, batch_rows = rows)), as.numeric(reference), tolerance = 1e-12)
  }
  expect_equal(as.numeric(el_vegdist(x, method = "aitchison", pseudocount = 0.5, threads = 4)), as.numeric(reference), tolerance = 1e-12)
})

test_that("a lazy stingy table is accepted for aitchison", {
  x <- duckplyr::as_duckdb_tibble(long_table(fixtures$synthetic), prudence = "stingy")

  expect_no_error(el_vegdist(x, method = "aitchison"))
})

test_that("bad aitchison input is an error", {
  x <- long_table(fixtures$synthetic)
  pc <- stats::setNames(rep(1, 3), c("s01", "s02", "s03"))

  expect_error(el_vegdist(x, method = "aitchison", pseudocount = 0), "`pseudocount` must be")
  expect_error(el_vegdist(x, method = "aitchison", pseudocount = -1), "`pseudocount` must be")
  expect_error(el_vegdist(x, method = "aitchison", pseudocount = "a"), "`pseudocount` must be")
  expect_error(el_vegdist(x, method = "aitchison", pseudocount = "relative"), "`pseudocount` must be")
  expect_error(el_vegdist(x, method = "aitchison", pseudocount = c(1, 2)), "`pseudocount` must be")
  expect_error(el_vegdist(x, method = "aitchison", pseudocount = pc), "no pseudocount for sample")
  expect_error(el_vegdist(transform(x, quant = -quant), method = "aitchison"), "negative")
  expect_error(el_vegdist(x, method = "aitchison", binary = TRUE), "binary")
})

for (fixture in c("varespec", "synthetic", "lmo")) {
  for (pc in c(1, 0.01)) {
    test_that(paste0("el_rda with transform = 'clr', pseudocount ", pc, " on ", fixture, " matches vegan::rda on decostand clr"), {
      m <- fixtures[[fixture]]
      got <- el_rda(long_table(m), k = 3, transform = "clr", pseudocount = pc, solver = "exact")

      compare_rda(got, vegan::decostand(m, "clr", pseudocount = pc))
    })
  }
}

test_that("el_rda clr with scale = TRUE matches vegan", {
  m <- fixtures$lmo
  got <- el_rda(long_table(m), k = 3, transform = "clr", scale = TRUE, solver = "exact")

  compare_rda(got, vegan::decostand(m, "clr", pseudocount = 1), scale = TRUE)
})

test_that("the iterative solver on clr agrees with vegan on a larger sparse table", {
  m <- present_taxa(big_sparse_community())
  got <- el_rda(long_table(m), k = 4, transform = "clr", pseudocount = 1, solver = "irlba")

  expect_equal(got$ecolongr$solver, "irlba")
  compare_rda(got, vegan::decostand(m, "clr", pseudocount = 1), k = 4, tolerance = 1e-6)
})

test_that("el_rda clr accepts a named pseudocount per sample and 'inverse_depth'", {
  m <- fixtures$lmo
  pc <- stats::setNames(1 / rowSums(m), rownames(m))

  for (value in list(pc, "inverse_depth")) {
    got <- el_rda(long_table(m), k = 3, transform = "clr", pseudocount = value, solver = "exact")
    compare_rda(got, clr_reference(m, pc))
  }
})

test_that("el_rda clr with an all-zero sample listed in `samples` matches vegan on the zero row", {
  m <- present_taxa(fixtures$synthetic[1:6, ])
  with_zero <- rbind(m, empty = 0)
  got <- el_rda(long_table(m), k = 3, transform = "clr", samples = rownames(with_zero), solver = "exact")

  compare_rda(got, vegan::decostand(with_zero, "clr", pseudocount = 1))
})

test_that("transform = 'none' is the default and ignores no argument silently", {
  m <- fixtures$varespec

  expect_equal(unname(vegan::eigenvals(el_rda(long_table(m), k = 3))), unname(vegan::eigenvals(el_rda(long_table(m), k = 3, transform = "none"))))
  expect_error(el_rda(long_table(m), pseudocount = 0.5), "`pseudocount` is only used with `transform = \"clr\"`")
  expect_error(el_rda(long_table(m), transform = "clr", pseudocount = 0), "`pseudocount` must be")
  expect_error(el_rda(long_table(m), transform = "hellinger"), "should be one of")
})

test_that("el_decostand points to the clr functions instead of building a dense table", {
  expect_error(el_decostand(long_table(fixtures$varespec), "clr"), "el_vegdist.*aitchison.*el_rda")
})

test_that("clr ordination of the LMO subset recovers the seasonal structure, as vegan's does", {
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
  got <- el_rda(lmo$long, k = 2, transform = "clr", pseudocount = 1, solver = "exact")
  reference <- vegan::rda(vegan::decostand(lmo$matrix, "clr", pseudocount = 1))

  expect_gt(seasonal_r2(got), 0.5)
  expect_equal(seasonal_r2(got), seasonal_r2(reference), tolerance = 1e-8)
})
