data(varespec, package = "vegan")

fixtures <- list(varespec = as.matrix(varespec), synthetic = synthetic_community())

for (fixture in names(fixtures)) {
  test_that(paste0("el_rda on ", fixture, " matches vegan::rda: eigenvalues, total inertia and scores for scalings 1, 2, 3"), {
    m <- fixtures[[fixture]]
    got <- el_rda(long_table(m), k = 3, solver = "exact")

    compare_rda(got, m)
  })
}

test_that("el_rda with scale = TRUE matches vegan::rda(scale = TRUE)", {
  m <- fixtures$varespec
  got <- el_rda(long_table(m), k = 3, scale = TRUE, solver = "exact")

  compare_rda(got, m, scale = TRUE)
})

test_that("the iterative solver agrees with vegan on a larger sparse table", {
  m <- big_sparse_community()
  got <- el_rda(long_table(m), k = 4, solver = "irlba")

  compare_rda(got, m, k = 4, tolerance = 1e-6)
})

test_that("the automatic solver picks the exact one for small and irlba for large tables", {
  small <- el_rda(long_table(fixtures$varespec), k = 2)
  large <- el_rda(long_table(big_sparse_community()), k = 2)

  expect_equal(small$ecolongr$solver, "exact")
  expect_equal(large$ecolongr$solver, "irlba")
})

test_that("both solvers give the same answer", {
  x <- long_table(big_sparse_community(n = 120, p = 400))
  a <- el_rda(x, k = 3, solver = "exact")
  b <- el_rda(x, k = 3, solver = "irlba")

  expect_equal(unname(vegan::eigenvals(a)), unname(vegan::eigenvals(b)), tolerance = 1e-6)
})

test_that("the site scores are orthonormal before scaling and the species loadings too", {
  got <- el_rda(long_table(fixtures$synthetic), k = 3, solver = "exact")

  expect_equal(crossprod(got$CA$u), diag(3), tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(crossprod(got$CA$v), diag(3), tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("an all-zero sample listed in `samples` counts in n, the means and the total inertia", {
  m <- fixtures$synthetic[1:6, ]
  with_zero <- rbind(m, empty = 0)
  got <- el_rda(long_table(m), k = 3, samples = rownames(with_zero), solver = "exact")

  compare_rda(got, with_zero)
  expect_equal(nrow(vegan::scores(got, display = "sites")), 7)
})

test_that("without `samples` the all-zero sample is not in the result", {
  m <- fixtures$synthetic[1:6, ]
  got <- el_rda(long_table(m), k = 3, solver = "exact")

  expect_equal(nrow(vegan::scores(got, display = "sites")), 6)
  expect_false("empty" %in% rownames(vegan::scores(got, display = "sites")))
})

test_that("a transformed lazy table gives vegan's result on the transformed matrix", {
  m <- fixtures$varespec

  got <- el_rda(el_decostand(long_table(m), "hellinger"), k = 3, solver = "exact")
  compare_rda(got, vegan::decostand(m, "hellinger"))

  got <- el_rda(el_decostand(long_table(m), "total"), k = 3, solver = "exact")
  compare_rda(got, vegan::decostand(m, "total"))
})

test_that("rclr with absent cells as 0 gives vegan's result on the rclr matrix with NA as 0", {
  for (m in list(fixtures$varespec, fixtures$synthetic)) {
    rclr <- vegan_decostand(m, "rclr")
    rclr[is.na(rclr)] <- 0
    got <- el_rda(el_decostand(long_table(m), "rclr"), k = 3, solver = "exact")

    compare_rda(got, rclr)
  }
})

test_that("negative quantities are accepted", {
  m <- fixtures$synthetic
  m[m > 0] <- m[m > 0] * sample(c(-1, 1), sum(m > 0), TRUE)
  got <- el_rda(long_table_nonzero(m), k = 2, solver = "exact")

  expect_gt(got$tot.chi, 0)
})

test_that("the result is a vegan rda object that vegan's tools accept", {
  got <- el_rda(long_table(fixtures$varespec), k = 3, solver = "exact")

  expect_s3_class(got, "rda")
  expect_s3_class(got, "cca")
  expect_equal(length(vegan::eigenvals(got)), 3)
  expect_equal(dim(vegan::scores(got, display = "sites")), c(24, 2))

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_no_error(plot(got))
  expect_no_error(stats::screeplot(got))
})

test_that("species scores are named by taxon and site scores by sample", {
  m <- fixtures$varespec
  got <- el_rda(long_table(m), k = 2, solver = "exact")

  expect_setequal(rownames(vegan::scores(got, display = "species")), colnames(m))
  expect_setequal(rownames(vegan::scores(got, display = "sites")), rownames(m))
})

test_that("k is capped at the rank and must be a positive whole number", {
  m <- fixtures$synthetic[1:5, ]
  got <- el_rda(long_table(m), k = 50, solver = "exact")

  expect_equal(length(vegan::eigenvals(got)), 4)
  expect_error(el_rda(long_table(m), k = 0), "`k` must be")
  expect_error(el_rda(long_table(m), k = 1.5), "`k` must be")
  expect_error(el_rda(long_table(m), k = "a"), "`k` must be")
})

test_that("a lazy stingy table is accepted", {
  x <- duckplyr::as_duckdb_tibble(long_table(fixtures$varespec), prudence = "stingy")

  expect_no_error(el_rda(x, k = 2))
})

test_that("bad input is an error", {
  x <- long_table(fixtures$synthetic)

  expect_error(el_rda(dplyr::rename(x, abundance = quant)), "quant")
  expect_error(el_rda(transform(x, quant = NA_real_)), "missing")
  expect_error(el_rda(x, solver = "magic"), "should be one of")
  constant <- data.frame(sample = rep(c("a", "b", "c"), 2), taxon = rep(c("t1", "t2"), each = 3), quant = 1)
  expect_error(el_rda(constant, scale = TRUE), "constant")
})

unexplained <- function(x) 1 - sum(unname(vegan::eigenvals(x))) / x$tot.chi

test_that("the unexplained variance with k below the rank equals vegan's, for none, clr and scale = TRUE", {
  m <- present_taxa(fixtures$synthetic)
  clr <- vegan::decostand(m, "clr", pseudocount = 1)
  cases <- list(
    list(got = el_rda(long_table(m), k = 3, solver = "exact"), ref = vegan::rda(m)),
    list(got = el_rda(long_table(m), k = 3, transform = "clr", pseudocount = 1, solver = "exact"), ref = vegan::rda(clr)),
    list(got = el_rda(long_table(m), k = 3, scale = TRUE, solver = "exact"), ref = vegan::rda(m, scale = TRUE))
  )

  for (case in cases) {
    ref_unexplained <- 1 - sum(vegan::eigenvals(case$ref)[1:3]) / case$ref$tot.chi
    expect_gt(ref_unexplained, 0.01)
    expect_equal(unexplained(case$got), ref_unexplained, tolerance = 1e-8)
  }
})

test_that("total inertia does not depend on k or on the solver, so only the unexplained share changes", {
  m <- present_taxa(big_sparse_community())
  few <- el_rda(long_table(m), k = 2, solver = "irlba")
  more <- el_rda(long_table(m), k = 8, solver = "irlba")
  exact <- el_rda(long_table(m), k = 8, solver = "exact")

  expect_equal(few$tot.chi, more$tot.chi)
  expect_equal(more$tot.chi, exact$tot.chi, tolerance = 1e-10)
  expect_equal(unexplained(more), unexplained(exact), tolerance = 1e-6)
  expect_lt(unexplained(more), unexplained(few))
  expect_equal(unexplained(few), 1 - sum(vegan::eigenvals(vegan::rda(m))[1:2]) / vegan::rda(m)$tot.chi, tolerance = 1e-6)
})

test_that("all axes explain everything: the unexplained share is zero at k = rank", {
  m <- fixtures$synthetic[1:6, ]
  got <- el_rda(long_table(m), k = 1000, solver = "exact")

  expect_equal(length(vegan::eigenvals(got)), 5)
  expect_equal(unexplained(got), 0, tolerance = 1e-10)
})

test_that("the three axes of the completed rclr explain all of its inertia", {
  m <- fixtures$varespec
  got <- el_rda(long_table(m), k = 3, transform = "rclr")

  expect_equal(unexplained(got), 0, tolerance = 1e-10)
})

test_that("the option ecolongr.quiet silences the report and unset it is off", {
  m <- present_taxa(fixtures$synthetic)

  withr::local_options(ecolongr.quiet = TRUE)
  expect_no_message(el_rda(long_table(m), k = 3, solver = "exact"), message = "left out")
  withr::local_options(ecolongr.quiet = NULL)
  expect_message(el_rda(long_table(m), k = 3, solver = "exact"), "left out")
})

test_that("el_rda reports the shares held and left out by the k axes, equal to vegan's", {
  withr::local_options(ecolongr.quiet = FALSE)
  m <- present_taxa(fixtures$synthetic)
  ref <- vegan::rda(m)
  held <- sum(vegan::eigenvals(ref)[1:3]) / ref$tot.chi
  expect_gt(1 - held, 0.01)

  expect_message(el_rda(long_table(m), k = 3, solver = "exact"), sprintf("%.1f%% of the total variance", 100 * held), fixed = TRUE)
  expect_message(el_rda(long_table(m), k = 3, solver = "exact"), sprintf("%.1f%% is left out", 100 * (1 - held)), fixed = TRUE)
})

test_that("the message names the number of axes returned, not the number asked for", {
  withr::local_options(ecolongr.quiet = FALSE)
  m <- fixtures$synthetic[1:5, ]

  expect_message(el_rda(long_table(m), k = 50, solver = "exact"), "first 4 axes")
})

test_that("the message says nothing is left out when the axes hold everything", {
  withr::local_options(ecolongr.quiet = FALSE)
  expect_message(el_rda(long_table(fixtures$varespec), k = 3, transform = "rclr"), "0.0% is left out", fixed = TRUE)
})

test_that("the message is a message, so suppressMessages silences it and the result is unchanged", {
  m <- present_taxa(fixtures$synthetic)
  quiet <- suppressMessages(el_rda(long_table(m), k = 3, solver = "exact"))

  expect_no_warning(suppressMessages(el_rda(long_table(m), k = 3, solver = "exact")))
  expect_equal(quiet$tot.chi, vegan::rda(m)$tot.chi, tolerance = 1e-10)
})

test_that("dividing the kept eigenvalues by their sum always gives 1, dividing by tot.chi gives the share of the whole", {
  m <- present_taxa(fixtures$synthetic)
  ref <- vegan::rda(m)
  got <- suppressMessages(el_rda(long_table(m), k = 3, solver = "exact"))
  eig <- unname(vegan::eigenvals(got))

  expect_equal(sum(eig / sum(eig)), 1)
  expect_lt(sum(eig / got$tot.chi), 0.99)
  expect_equal(as.numeric(eig / got$tot.chi), unname(vegan::eigenvals(ref)[1:3] / ref$tot.chi), tolerance = 1e-8)
})

test_that("summary() gives vegan's proportion of the whole for the kept axes, not of their sum", {
  m <- present_taxa(fixtures$synthetic)
  ref_importance <- summary(vegan::rda(m))$cont$importance
  got <- suppressMessages(el_rda(long_table(m), k = 3, solver = "exact"))
  got_importance <- summary(got)$cont$importance

  expect_equal(unname(got_importance["Proportion Explained", 1:3]), unname(ref_importance["Proportion Explained", 1:3]), tolerance = 1e-8)
  expect_lt(sum(got_importance["Proportion Explained", 1:3]), 0.99)
})
