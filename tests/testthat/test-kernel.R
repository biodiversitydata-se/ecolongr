naive_pairs <- function(sample, taxon, quant, n, product = FALSE) {
  acc <- numeric(n * (n - 1) / 2)
  for (t in unique(taxon)) {
    rows <- which(taxon == t)
    for (a in rows) for (b in rows) {
      if (sample[a] < sample[b]) {
        k <- sample[a] * (2 * n - sample[a] - 1) / 2 + (sample[b] - sample[a] - 1) + 1
        acc[k] <- acc[k] + if (product) quant[a] * quant[b] else min(quant[a], quant[b])
      }
    }
  }
  acc
}

pairs_table <- function() {
  data.frame(
    sample = c(0L, 2L, 3L, 1L, 2L, 0L, 1L, 2L, 3L, 4L, 4L),
    taxon  = c(1L, 1L, 1L, 2L, 2L, 3L, 3L, 3L, 3L, 3L, 4L),
    quant  = c(5, 2, 9, 4, 1, 3, 3, 8, 6, 2, 7)
  )
}

test_that("accumulate_pairs adds the minimum over shared taxa in dist order", {
  d <- pairs_table()
  acc <- numeric(5 * 4 / 2)
  accumulate_pairs(acc, d$sample, d$taxon, d$quant, 5L, product = FALSE)

  expect_equal(acc, naive_pairs(d$sample, d$taxon, d$quant, 5L))
})

test_that("accumulate_pairs adds the product when asked", {
  d <- pairs_table()
  acc <- numeric(5 * 4 / 2)
  accumulate_pairs(acc, d$sample, d$taxon, d$quant, 5L, product = TRUE)

  expect_equal(acc, naive_pairs(d$sample, d$taxon, d$quant, 5L, product = TRUE))
})

test_that("accumulate_pairs adds to what is already in the accumulator", {
  d <- pairs_table()
  first <- d$taxon <= 2
  acc <- numeric(10)
  accumulate_pairs(acc, d$sample[first], d$taxon[first], d$quant[first], 5L, product = FALSE)
  accumulate_pairs(acc, d$sample[!first], d$taxon[!first], d$quant[!first], 5L, product = FALSE)

  expect_equal(acc, naive_pairs(d$sample, d$taxon, d$quant, 5L))
})

test_that("accumulate_pairs modifies the accumulator in place", {
  d <- pairs_table()
  acc <- numeric(10)
  alias <- acc
  accumulate_pairs(acc, d$sample, d$taxon, d$quant, 5L, product = FALSE)

  expect_equal(sum(alias), sum(acc))
  expect_gt(sum(acc), 0)
})

test_that("a taxon in a single sample adds nothing", {
  acc <- numeric(3)
  accumulate_pairs(acc, c(0L, 1L, 2L), c(1L, 2L, 3L), c(4, 5, 6), 3L, product = FALSE)

  expect_equal(acc, c(0, 0, 0))
})

shuffled_within_taxon <- function(d, seed = 3) {
  set.seed(seed)
  d[order(d$taxon, runif(nrow(d))), ]
}

test_that("accumulate_pairs does not need the samples sorted within a taxon", {
  d <- shuffled_within_taxon(pairs_table())
  expect_false(identical(d$sample, pairs_table()$sample))

  for (product in c(FALSE, TRUE)) {
    acc <- numeric(10)
    accumulate_pairs(acc, d$sample, d$taxon, d$quant, 5L, product = product)
    expect_equal(acc, naive_pairs(d$sample, d$taxon, d$quant, 5L, product = product))
  }
})

random_community <- function(n_samples = 60L, n_taxa = 300L, seed = 11) {
  set.seed(seed)
  d <- unique(data.frame(
    sample = sample.int(n_samples, 1500, TRUE) - 1L,
    taxon = sort(sample.int(n_taxa, 1500, TRUE))
  ))
  d$quant <- rpois(nrow(d), 5) + 1
  shuffled_within_taxon(d)
}

test_that("accumulate_pairs gives identical results for any number of threads", {
  d <- random_community()
  n <- 60L
  reference <- numeric(n * (n - 1) / 2)
  accumulate_pairs(reference, d$sample, d$taxon, d$quant, n, product = FALSE, threads = 1L)

  expect_gt(sum(reference), 0)
  for (threads in c(2L, 3L, 7L, 64L)) {
    acc <- numeric(n * (n - 1) / 2)
    accumulate_pairs(acc, d$sample, d$taxon, d$quant, n, product = FALSE, threads = threads)
    expect_identical(acc, reference)
  }
})

test_that("accumulate_pairs with threads matches the reference loop", {
  d <- random_community(n_samples = 30L, n_taxa = 80L, seed = 5)
  acc <- numeric(30 * 29 / 2)
  accumulate_pairs(acc, d$sample, d$taxon, d$quant, 30L, product = TRUE, threads = 4L)

  expect_equal(acc, naive_pairs(d$sample, d$taxon, d$quant, 30L, product = TRUE))
})

test_that("accumulate_pairs with more threads than samples still works", {
  d <- pairs_table()
  acc <- numeric(10)
  accumulate_pairs(acc, d$sample, d$taxon, d$quant, 5L, product = FALSE, threads = 16L)

  expect_equal(acc, naive_pairs(d$sample, d$taxon, d$quant, 5L))
})

test_that("a table spanning many cache blocks matches vegan", {
  set.seed(21)
  m <- matrix(0, 1500, 300, dimnames = list(sprintf("s%04d", 1:1500), sprintf("t%03d", 1:300)))
  m[sample(length(m), 40000)] <- rpois(40000, 3) + 1
  m <- m[rowSums(m) > 0, ]

  for (threads in c(1L, 5L)) {
    got <- el_vegdist(long_table(m), threads = threads)
    expect_equal(as_full_matrix(got), as_full_matrix(vegan::vegdist(m)), tolerance = 1e-10)
  }
})
