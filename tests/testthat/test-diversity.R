data(BCI, package = "vegan")
data(varespec, package = "vegan")

diversity_fixtures <- list(
  varespec = as.matrix(varespec),
  BCI = as.matrix(BCI[1:15, ]),
  synthetic = synthetic_community()
)

vegan_diversity <- function(m, index, base = exp(1)) c(vegan::diversity(m, index, base = base))

for (fixture in names(diversity_fixtures)) {
  m <- diversity_fixtures[[fixture]]

  for (index in c("shannon", "simpson", "invsimpson")) {
    test_that(paste0("el_diversity(", index, ") matches vegan::diversity on ", fixture), {
      got <- el_diversity(long_table(m), index = index)

      expect_equal(got[rownames(m)], vegan_diversity(m, index), tolerance = 1e-10)
    })
  }

  for (base in c(2, 10)) {
    test_that(paste0("el_diversity(shannon) with base ", base, " matches vegan on ", fixture), {
      got <- el_diversity(long_table(m), "shannon", base = base)

      expect_equal(got[rownames(m)], vegan_diversity(m, "shannon", base), tolerance = 1e-10)
    })
  }

  test_that(paste0("el_specnumber matches vegan::specnumber on ", fixture), {
    got <- el_specnumber(long_table(m))

    expect_equal(got[rownames(m)], c(vegan::specnumber(m)), tolerance = 1e-10)
  })
}

test_that("the result is a numeric vector named by sample, one entry per sample", {
  m <- diversity_fixtures$synthetic
  for (got in list(el_diversity(long_table(m)), el_specnumber(long_table(m)))) {
    expect_true(is.numeric(got))
    expect_null(dim(got))
    expect_setequal(names(got), rownames(m))
    expect_length(got, nrow(m))
  }
})

test_that("a sample with one taxon has Shannon 0, Simpson 0, inverse Simpson 1 and richness 1", {
  x <- data.frame(sample = c("a", "b", "b"), taxon = c("t1", "t1", "t2"), quant = c(7, 3, 3))

  expect_equal(el_diversity(x, "shannon")[["a"]], 0)
  expect_equal(el_diversity(x, "simpson")[["a"]], 0)
  expect_equal(el_diversity(x, "invsimpson")[["a"]], 1)
  expect_equal(el_specnumber(x)[["a"]], 1)
  expect_equal(el_diversity(x, "shannon")[["b"]], log(2))
  expect_equal(el_diversity(x, "simpson")[["b"]], 0.5)
  expect_equal(el_diversity(x, "invsimpson")[["b"]], 2)
})

test_that("identical samples get identical values and are not merged", {
  m <- diversity_fixtures$synthetic
  expect_equal(m[1, ], m[2, ])

  for (got in list(el_diversity(long_table(m)), el_specnumber(long_table(m)))) {
    expect_equal(got[["s01"]], got[["s02"]])
  }
  expect_equal(el_specnumber(long_table(m))[["s01"]], sum(m[1, ] > 0))
})

test_that("proportions are accepted and a common scale does not change the result", {
  m <- diversity_fixtures$varespec
  p <- m / rowSums(m)
  scaled <- long_table(m)
  scaled$quant <- scaled$quant * 1000

  for (index in c("shannon", "simpson", "invsimpson")) {
    expect_equal(el_diversity(long_table(p), index)[rownames(m)], vegan_diversity(p, index), tolerance = 1e-10)
    expect_equal(el_diversity(scaled, index)[rownames(m)], vegan_diversity(m, index), tolerance = 1e-10)
  }
  expect_equal(el_specnumber(long_table(p))[rownames(m)], c(vegan::specnumber(p)), tolerance = 1e-10)
})

test_that("fractional quantities that are not proportions are normalised per sample", {
  x <- data.frame(sample = c("a", "a", "b", "b"), taxon = c("t1", "t2", "t1", "t2"), quant = c(0.25, 0.25, 30, 10))

  expect_equal(el_diversity(x, "simpson")[["a"]], 0.5)
  expect_equal(el_diversity(x, "simpson")[["b"]], 1 - (0.75^2 + 0.25^2))
})

test_that("the order of the rows does not matter", {
  x <- long_table(diversity_fixtures$BCI)
  set.seed(3)
  shuffled <- x[sample(nrow(x)), ]

  expect_equal(el_diversity(shuffled)[names(el_diversity(x))], el_diversity(x), tolerance = 1e-10)
  expect_equal(el_specnumber(shuffled)[names(el_specnumber(x))], el_specnumber(x))
})

test_that("rows with quant == 0 are ignored by both functions", {
  m <- diversity_fixtures$synthetic
  x <- long_table(m)
  absent <- which(m == 0, arr.ind = TRUE)[1:40, ]
  zero <- data.frame(sample = rownames(m)[absent[, 1]], taxon = colnames(m)[absent[, 2]], quant = 0)
  with_zero <- rbind(x, zero)

  expect_equal(el_specnumber(with_zero)[rownames(m)], c(vegan::specnumber(m)))
  for (index in c("shannon", "simpson", "invsimpson")) {
    expect_equal(el_diversity(with_zero, index)[rownames(m)], vegan_diversity(m, index), tolerance = 1e-10)
  }
})

test_that("a sample that has only zero rows is absent, and vegan gives it Shannon 0, Simpson 1, inverse Simpson Inf and richness 0", {
  m <- rbind(diversity_fixtures$varespec[1:3, ], empty = 0)
  x <- rbind(long_table(m), data.frame(sample = "empty", taxon = colnames(m)[1:2], quant = 0))

  expect_false("empty" %in% names(el_diversity(x)))
  expect_false("empty" %in% names(el_specnumber(x)))
  expect_length(el_specnumber(x), 3)

  expect_equal(vegan_diversity(m, "shannon")[["empty"]], 0)
  expect_equal(vegan_diversity(m, "simpson")[["empty"]], 1)
  expect_equal(vegan_diversity(m, "invsimpson")[["empty"]], Inf)
  expect_equal(c(vegan::specnumber(m))[["empty"]], 0)
})

test_that("negative and missing quantities are errors that name the column", {
  x <- long_table(diversity_fixtures$varespec[1:3, ])
  names(x)[3] <- "n"
  for (f in list(el_diversity, el_specnumber)) {
    y <- x
    y$n[1] <- -1
    expect_error(f(y, quant = "n"), "`n` has negative values")
    y$n[1] <- NA
    expect_error(f(y, quant = "n"), "`n` has missing values")
  }
})

test_that("check = FALSE skips the quantity check", {
  x <- long_table(diversity_fixtures$varespec[1:3, ])
  x$quant[1] <- -1

  expect_no_error(el_diversity(x, check = FALSE))
  expect_no_error(el_specnumber(x, check = FALSE))
})

test_that("an unknown index and a bad base are errors", {
  x <- long_table(diversity_fixtures$varespec[1:3, ])

  expect_error(el_diversity(x, "chao"), "should be one of")
  for (base in list(1, 0, -2, NA, Inf, c(2, 10), "2")) {
    expect_error(el_diversity(x, "shannon", base = base), "base")
  }
})

test_that("column names are taken from the arguments", {
  x <- long_table(diversity_fixtures$varespec[1:3, ])
  expected <- el_diversity(x)
  expected_richness <- el_specnumber(x)
  names(x) <- c("site", "species", "n")

  expect_equal(el_diversity(x, sample = "site", taxon = "species", quant = "n"), expected)
  expect_equal(el_specnumber(x, sample = "site", taxon = "species", quant = "n"), expected_richness)
})

test_that("both functions convert their input with prudence 'stingy'", {
  seen <- NULL
  original <- duckplyr::as_duckdb_tibble
  testthat::local_mocked_bindings(
    as_duckdb_tibble = function(x, ..., prudence = "lavish") {
      seen <<- c(seen, prudence)
      original(x, ..., prudence = prudence)
    },
    .package = "duckplyr"
  )
  x <- long_table(diversity_fixtures$varespec[1:3, ])

  el_diversity(x)
  el_specnumber(x)

  expect_gte(length(seen), 1)
  expect_true(all(seen == "stingy"))
})

test_that("base changes only Shannon", {
  x <- long_table(diversity_fixtures$varespec)

  expect_equal(el_diversity(x, "simpson", base = 2), el_diversity(x, "simpson"))
  expect_equal(el_diversity(x, "invsimpson", base = 10), el_diversity(x, "invsimpson"))
  expect_false(isTRUE(all.equal(el_diversity(x, "shannon", base = 2), el_diversity(x, "shannon"))))
})

test_that("the result is sorted by sample name whatever the order of the rows", {
  x <- data.frame(sample = c("b", "c", "a", "b", "c", "a"), taxon = rep(c("t1", "t2"), each = 3), quant = c(1, 2, 3, 4, 5, 6))

  expect_identical(names(el_diversity(x)), c("a", "b", "c"))
  expect_identical(names(el_specnumber(x)), c("a", "b", "c"))
  expect_equal(el_specnumber(x), c(a = 2, b = 2, c = 2))
  expect_equal(el_diversity(x, "simpson")[["a"]], 1 - (3 / 9)^2 - (6 / 9)^2)
  expect_equal(el_diversity(x, "simpson")[["c"]], 1 - (2 / 7)^2 - (5 / 7)^2)
})

test_that("sample names are the samples as character: factors sort as text, an NA sample is named NA, numbers keep their digits", {
  f <- data.frame(sample = factor(c("b", "a", "c"), levels = c("c", "b", "a")), taxon = "t1", quant = 1)
  expect_identical(names(el_specnumber(f)), c("a", "b", "c"))

  n <- data.frame(sample = c(10, 9, NA), taxon = "t1", quant = 1)
  got <- el_specnumber(n)
  expect_identical(names(got), c("10", "9", NA))
  expect_equal(unname(got), c(1, 1, 1))

  n$taxon <- c("t1", "t2", "t1")
  n$quant <- c(1, 1, 4)
  expect_identical(names(el_diversity(n, "simpson")), c("10", "9", NA))
  expect_equal(unname(el_diversity(n, "simpson")), c(0, 0, 0))
})

test_that("each row is a cell: a duplicated sample and taxon is not merged", {
  x <- data.frame(sample = c("a", "a"), taxon = c("t1", "t1"), quant = c(5, 5))

  expect_equal(el_specnumber(x)[["a"]], 2)
  expect_equal(el_diversity(x, "simpson")[["a"]], 0.5)
})

test_that("NaN and Inf quantities are errors that name the column", {
  x <- long_table(diversity_fixtures$varespec[1:3, ])
  names(x)[3] <- "n"
  for (bad in c(NaN, Inf)) {
    y <- x
    y$n[1] <- bad
    expect_error(el_diversity(y, quant = "n"), "`n` has non-finite values")
    expect_error(el_specnumber(y, quant = "n"), "`n` has non-finite values")
  }
})
