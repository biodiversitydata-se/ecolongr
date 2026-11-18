data(BCI, package = "vegan")

rarefy_fixtures <- list(BCI = as.matrix(BCI[1:12, ]), synthetic = synthetic_community()[-3, ])

vegan_drarefy <- function(m, depth) {
  out <- suppressWarnings(vegan::drarefy(m, depth))
  attributes(out) <- attributes(m)
  out
}

vegan_rarefy <- function(m, depth) {
  out <- suppressWarnings(vegan::rarefy(m, depth))
  c(out)
}

for (fixture in names(rarefy_fixtures)) {
  m <- rarefy_fixtures[[fixture]]
  m <- m[rowSums(m) > 0, colSums(m) > 0]
  smallest <- min(rowSums(m))

  for (depth in unique(c(1, 5, floor(smallest / 2), smallest - 1, smallest))) {
    test_that(paste0("el_drarefy matches vegan::drarefy on ", fixture, " at depth ", depth), {
      got <- as.data.frame(el_drarefy(long_table(m), depth = depth))

      expect_equal(wide_from_long(got, m), vegan_drarefy(m, depth), tolerance = 1e-10)
      expect_equal(nrow(got), sum(m > 0))
    })

    test_that(paste0("el_rarefy matches vegan::rarefy on ", fixture, " at depth ", depth), {
      got <- el_rarefy(long_table(m), depth = depth)

      expect_equal(got[rownames(m)], vegan_rarefy(m, depth), tolerance = 1e-10)
    })
  }
}

test_that("a depth above a sample total gives probability 1 and the observed richness, as vegan", {
  m <- rarefy_fixtures$BCI[1:4, ]
  depth <- max(rowSums(m)) + 10

  expect_warning(got <- as.data.frame(el_drarefy(long_table(m), depth = depth)), "not rarefied")
  expect_equal(got$quant, rep(1, nrow(got)))
  expect_equal(wide_from_long(got, m), vegan_drarefy(m, depth), tolerance = 1e-10)
  expect_warning(richness <- el_rarefy(long_table(m), depth = depth), "not rarefied")
  expect_equal(richness[rownames(m)], rowSums(m > 0), ignore_attr = TRUE)
  expect_equal(richness[rownames(m)], vegan_rarefy(m, depth), tolerance = 1e-10)
})

test_that("a depth equal to every total gives no warning", {
  m <- matrix(c(2, 3, 5, 4, 4, 2), 2, byrow = TRUE, dimnames = list(c("a", "b"), c("x", "y", "z")))

  expect_warning(el_rarefy(long_table(m), depth = 10), NA)
})

test_that("the expected values are the limit of el_rrarefy over many replicates", {
  m <- rarefy_fixture()
  depth <- 11
  reps <- 4000
  set.seed(5)
  drawn <- as.data.frame(el_rrarefy(replicated(m, reps), depth = depth))
  drawn$base <- sub("_.*$", "", drawn$sample)

  observed_richness <- tapply(drawn$taxon, drawn$sample, length)
  mean_richness <- tapply(observed_richness, sub("_.*$", "", names(observed_richness)), mean)
  sd_richness <- tapply(observed_richness, sub("_.*$", "", names(observed_richness)), stats::sd)
  expected <- el_rarefy(long_table(m), depth = depth)
  expect_true(all(abs(mean_richness[names(expected)] - expected) < 5 * sd_richness[names(expected)] / sqrt(reps)))

  present <- tapply(rep(1, nrow(drawn)), list(drawn$base, drawn$taxon), sum) / reps
  probability <- as.data.frame(el_drarefy(long_table(m), depth = depth))
  p <- probability$quant
  got <- present[cbind(probability$sample, probability$taxon)]
  expect_true(all(abs(got - p) < 5 * sqrt(p * (1 - p) / reps) + 1e-12))
})

test_that("non-integer, negative and missing counts and a bad depth are errors", {
  x <- long_table(rarefy_fixture())
  for (f in list(el_drarefy, el_rarefy)) {
    y <- x
    y$quant[1] <- 1.5
    expect_error(f(y, depth = 5), "integer")
    y$quant[1] <- -1
    expect_error(f(y, depth = 5), "negative")
    y$quant[1] <- NA
    expect_error(f(y, depth = 5), "missing")
    expect_error(f(x, depth = 0), "depth")
    expect_error(f(x, depth = c(5, 6)), "depth")
    expect_error(f(x, depth = Inf), "depth")
  }
})

test_that("column names are taken from the arguments", {
  x <- long_table(rarefy_fixture())
  names(x) <- c("site", "species", "n")
  got <- as.data.frame(el_drarefy(x, depth = 5, sample = "site", taxon = "species", quant = "n"))

  expect_named(got, c("site", "species", "n"))
  expect_setequal(names(el_rarefy(x, depth = 5, sample = "site", taxon = "species", quant = "n")), c("a", "b", "c"))
})

test_that("rows with a zero count are ignored", {
  m <- rarefy_fixture()
  x <- long_table(m)
  zero <- data.frame(sample = "a", taxon = "t6", quant = 0)

  sorted <- function(d) d[order(d$sample, d$taxon), ]
  with_zero <- sorted(as.data.frame(el_drarefy(rbind(x, zero), depth = 5)))
  without <- sorted(as.data.frame(el_drarefy(x, depth = 5)))

  expect_equal(with_zero, without, ignore_attr = TRUE)
  expect_equal(nrow(with_zero), nrow(x))
})

test_that("el_drarefy returns a lazy duckplyr table and both convert with prudence 'stingy'", {
  seen <- NULL
  original <- duckplyr::as_duckdb_tibble
  testthat::local_mocked_bindings(
    as_duckdb_tibble = function(x, ..., prudence = "lavish") {
      seen <<- c(seen, prudence)
      original(x, ..., prudence = prudence)
    },
    .package = "duckplyr"
  )
  x <- long_table(rarefy_fixture())

  expect_s3_class(el_drarefy(x, depth = 5), "duckplyr_df")
  el_rarefy(x, depth = 5)

  expect_true(length(seen) > 0)
  expect_true(all(seen == "stingy"))
})

test_that("large totals and a small depth still match vegan to 1e-10", {
  set.seed(11)
  counts <- as.numeric(rmultinom(1, 2e6, prob = 1 / (1:500)))
  m <- matrix(counts, 1, dimnames = list("a", paste0("t", 1:500)))
  x <- long_table(m)

  for (depth in c(10, 1000, 1e5)) {
    expect_equal(el_rarefy(x, depth)[["a"]], vegan_rarefy(m, depth)[[1]], tolerance = 1e-10)
    got <- as.data.frame(el_drarefy(x, depth))
    expect_equal(got$quant, as.numeric(vegan_drarefy(m, depth))[match(got$taxon, colnames(m))], tolerance = 1e-10)
  }
})
