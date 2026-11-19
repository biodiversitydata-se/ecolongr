cell_moments <- function(got, m, reps) {
  key <- rep(rownames(m), each = ncol(m))
  full <- matrix(0, nrow(m) * reps, ncol(m), dimnames = list(NULL, colnames(m)))
  parts <- strsplit(got$sample, "_", fixed = TRUE)
  row <- (as.integer(vapply(parts, `[`, "", 2)) - 1L) * nrow(m) + match(vapply(parts, `[`, "", 1), rownames(m))
  full[cbind(row, match(got$taxon, colnames(m)))] <- got$quant
  sample_of_row <- rep(rownames(m), times = reps)
  list(
    mean = t(sapply(rownames(m), function(s) colMeans(full[sample_of_row == s, , drop = FALSE]))),
    var = t(sapply(rownames(m), function(s) apply(full[sample_of_row == s, , drop = FALSE], 2, stats::var)))
  )
}

hypergeometric_moments <- function(m, depth) {
  N <- rowSums(m)
  p <- m / N
  list(mean = depth * p, var = depth * p * (1 - p) * (N - depth) / (N - 1))
}

test_that("every sample is rarefied to exactly the depth, using only counts it had", {
  m <- rarefy_fixture()
  set.seed(1)
  got <- as.data.frame(el_rrarefy(long_table(m), depth = 11))

  expect_equal(c(tapply(got$quant, got$sample, sum)), c(a = 11, b = 11, c = 11))
  expect_true(all(got$quant >= 1))
  expect_equal(got$quant, round(got$quant))
  expect_true(all(got$quant <= m[cbind(match(got$sample, rownames(m)), match(got$taxon, colnames(m)))]))
  expect_false(any(got$taxon == "t6"))
})

test_that("rarefied counts have the hypergeometric mean and variance, as in vegan::rrarefy", {
  m <- rarefy_fixture()
  depth <- 11
  reps <- 4000
  exact <- hypergeometric_moments(m, depth)

  set.seed(2)
  ours <- cell_moments(as.data.frame(el_rrarefy(replicated(m, reps), depth = depth)), m, reps)
  set.seed(3)
  theirs_matrix <- vegan::rrarefy(m[rep(rownames(m), times = reps), ], depth)
  sample_of_row <- rep(rownames(m), times = reps)
  theirs <- list(
    mean = t(sapply(rownames(m), function(s) colMeans(theirs_matrix[sample_of_row == s, , drop = FALSE]))),
    var = t(sapply(rownames(m), function(s) apply(theirs_matrix[sample_of_row == s, , drop = FALSE], 2, stats::var)))
  )

  present <- m > 0
  se <- sqrt(exact$var / reps)
  expect_true(all(abs(ours$mean - exact$mean)[present] < 5 * se[present]))
  expect_true(all(abs(theirs$mean - exact$mean)[present] < 5 * se[present]))
  expect_true(all(abs(ours$mean - theirs$mean)[present] < 5 * sqrt(2) * se[present]))
  expect_true(all(abs(ours$var / exact$var - 1)[present] < 0.2))
  expect_true(all(abs(theirs$var / exact$var - 1)[present] < 0.2))
  expect_true(all(abs(ours$var / theirs$var - 1)[present] < 0.2))
})

test_that("a sample whose total equals the depth is kept as it is, without a warning", {
  m <- rarefy_fixture() # totals 20, 30, 12

  expect_warning(got <- as.data.frame(el_rrarefy(long_table(m), depth = 12)), NA)
  expect_equal(got[got$sample == "c", c("taxon", "quant")][order(got$taxon[got$sample == "c"]), ],
               long_table(m)[long_table(m)$sample == "c", c("taxon", "quant")],
               ignore_attr = TRUE)
  expect_equal(c(tapply(got$quant, got$sample, sum)), c(a = 12, b = 12, c = 12))
})

test_that("samples with a total below the depth are dropped, with a warning that names them", {
  m <- rarefy_fixture() # totals 20, 30, 12
  m <- rbind(m, d = c(1, 1, 0, 0, 0, 0))

  expect_warning(got <- as.data.frame(el_rrarefy(long_table(m), depth = 15)), "dropped.*c, d")
  expect_equal(sort(unique(got$sample)), c("a", "b"))
  expect_equal(c(tapply(got$quant, got$sample, sum)), c(a = 15, b = 15))
})

test_that("the warning counts the dropped samples and shortens a long list of names", {
  m <- matrix(1, nrow = 12, ncol = 2, dimnames = list(sprintf("s%02d", 1:12), c("x", "y")))
  m <- rbind(m, big = c(50, 50))

  expect_warning(el_rrarefy(long_table(m), depth = 10), "12 samples.*s01.*\\.\\.\\.")
})

test_that("it is an error when no sample reaches the depth", {
  expect_error(el_rrarefy(long_table(rarefy_fixture()), depth = 31), "depth")
})

test_that("without a depth, samples are rarefied to the smallest total, with no warning", {
  m <- rarefy_fixture() # totals 20, 30, 12
  set.seed(4)

  expect_warning(got <- as.data.frame(el_rrarefy(long_table(m))), NA)
  expect_equal(c(tapply(got$quant, got$sample, sum)), c(a = 12, b = 12, c = 12))
})

test_that("the smallest total ignores zero counts and samples that have only zeros", {
  x <- data.frame(sample = c("a", "a", "b", "b", "z"), taxon = c("x", "y", "x", "y", "x"),
                  quant = c(5, 0, 8, 9, 0))

  expect_warning(got <- as.data.frame(el_rrarefy(x)), NA)
  expect_equal(c(tapply(got$quant, got$sample, sum)), c(a = 5, b = 5))
})

test_that("the same seed gives the same result and another seed another", {
  x <- long_table(rarefy_fixture())
  set.seed(7)
  one <- as.data.frame(el_rrarefy(x, depth = 11))
  set.seed(7)
  two <- as.data.frame(el_rrarefy(x, depth = 11))
  set.seed(8)
  three <- as.data.frame(el_rrarefy(x, depth = 11))

  expect_equal(one, two)
  expect_false(isTRUE(all.equal(one, three)))
})

test_that("non-integer counts and a bad depth are errors", {
  x <- long_table(rarefy_fixture())
  y <- x
  y$quant[1] <- 1.5

  expect_error(el_rrarefy(y, depth = 5), "integer")
  expect_error(el_rrarefy(x, depth = 0), "depth")
  expect_error(el_rrarefy(x, depth = c(5, 6)), "depth")
  expect_error(el_rrarefy(x, depth = NA), "depth")
  expect_error(el_rrarefy(x, depth = 5.5), "depth")
})

test_that("column names are taken from the arguments and kept in the result", {
  x <- long_table(rarefy_fixture())
  names(x) <- c("site", "species", "n")
  got <- as.data.frame(el_rrarefy(x, depth = 11, sample = "site", taxon = "species", quant = "n"))

  expect_named(got, c("site", "species", "n"))
  expect_equal(sum(got$n[got$site == "a"]), 11)
})

test_that("a lazy duckplyr table is accepted and the result is a duckplyr table", {
  x <- duckplyr::as_duckdb_tibble(long_table(rarefy_fixture()), prudence = "stingy")
  got <- el_rrarefy(x, depth = 11)

  expect_s3_class(got, "duckplyr_df")
  expect_equal(sum(as.data.frame(got)$quant), 33)
})

test_that("missing, negative, infinite and NaN counts are errors that name the column", {
  x <- long_table(rarefy_fixture())
  names(x)[3] <- "n"
  for (bad in list(NA_real_, -1, Inf, NaN)) {
    y <- x
    y$n[1] <- bad
    expect_error(el_rrarefy(y, depth = 5, quant = "n"), "`n`")
  }
})

test_that("depth must be finite", {
  expect_error(el_rrarefy(long_table(rarefy_fixture()), depth = Inf), "depth")
})

test_that("counts whose total exceeds the integer range are rarefied", {
  x <- data.frame(sample = "a", taxon = c("x", "y"), quant = c(2000000000L, 2000000000L))
  got <- as.data.frame(el_rrarefy(x, depth = 5))

  expect_equal(sum(got$quant), 5)
})

test_that("a sample is one group even when other ids collate as equal and rows are interleaved", {
  nfc <- "á"
  nfd <- "á"
  x <- data.frame(sample = c(nfc, nfd, nfc), taxon = c("x", "x", "y"), quant = c(5, 5, 5))
  expect_warning(got <- as.data.frame(el_rrarefy(x, depth = 6)), "dropped")

  expect_equal(sum(got$quant[got$sample == nfc]), 6)
  expect_false(nfd %in% got$sample)
})

test_that("rows with a zero count are ignored, and a sample with only zeros is dropped without a warning", {
  x <- data.frame(sample = c("a", "a", "b"), taxon = c("x", "y", "x"), quant = c(5, 0, 0))

  expect_warning(got <- as.data.frame(el_rrarefy(x, depth = 3)), NA)
  expect_equal(got, data.frame(sample = "a", taxon = "x", quant = 3))
})

test_that("samples come out in the order they first appear", {
  x <- data.frame(sample = c("b", "a", "b"), taxon = c("x", "x", "y"), quant = c(5, 5, 5))
  got <- as.data.frame(el_rrarefy(x, depth = 3))

  expect_equal(unique(got$sample), c("b", "a"))
})

test_that("el_rrarefy converts its input with prudence 'stingy'", {
  seen <- NULL
  original <- duckplyr::as_duckdb_tibble
  testthat::local_mocked_bindings(
    as_duckdb_tibble = function(x, ..., prudence = "lavish") {
      seen <<- c(seen, prudence)
      original(x, ..., prudence = prudence)
    },
    .package = "duckplyr"
  )

  el_rrarefy(long_table(rarefy_fixture()), depth = 5)

  expect_true(length(seen) > 0)
  expect_true(all(seen == "stingy"))
})
