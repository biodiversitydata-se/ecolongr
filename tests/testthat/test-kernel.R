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
