data(BCI, package = "vegan")

frequency_fixtures <- list(BCI = as.matrix(BCI[1:15, ]), synthetic = synthetic_community())

reference_frequencies <- function(m) {
  do.call(rbind, lapply(rownames(m), function(s) {
    counts <- m[s, m[s, ] > 0]
    if (length(counts) == 0) return(NULL)
    tab <- table(counts)
    data.frame(sample = s, count = as.numeric(names(tab)), n_taxa = as.numeric(tab))
  }))
}

sorted_frequencies <- function(d) {
  d <- data.frame(sample = d[[1]], count = as.numeric(d$count), n_taxa = as.numeric(d$n_taxa))
  d <- d[order(d$sample, d$count), ]
  rownames(d) <- NULL
  d
}

for (fixture in names(frequency_fixtures)) {
  test_that(paste0("el_frequency_counts matches a table of the counts of each sample on ", fixture), {
    m <- frequency_fixtures[[fixture]]
    got <- sorted_frequencies(as.data.frame(el_frequency_counts(long_table(m))))

    expect_equal(got, sorted_frequencies(reference_frequencies(m)))
  })

  test_that(paste0("the frequencies add up to the richness and the total of each sample on ", fixture), {
    m <- frequency_fixtures[[fixture]]
    m <- m[rowSums(m) > 0, ]
    got <- as.data.frame(el_frequency_counts(long_table(m)))

    expect_equal(c(tapply(got$n_taxa, got$sample, sum))[rownames(m)], vegan::specnumber(m), ignore_attr = TRUE)
    expect_equal(c(tapply(got$count * got$n_taxa, got$sample, sum))[rownames(m)], rowSums(m), ignore_attr = TRUE)
  })
}

test_that("a small table gives the hand-counted frequencies, with samples kept apart", {
  x <- data.frame(
    sample = c("a", "a", "a", "a", "b", "b", "b"),
    taxon  = c("t1", "t2", "t3", "t4", "t1", "t2", "t3"),
    quant  = c(1, 1, 2, 5, 2, 2, 2)
  )
  got <- sorted_frequencies(as.data.frame(el_frequency_counts(x)))

  expect_equal(got, data.frame(sample = c("a", "a", "a", "b"), count = c(1, 2, 5, 2), n_taxa = c(2, 1, 1, 3)))
})

test_that("rows with a zero count are ignored and a sample with only zeros is absent", {
  x <- data.frame(sample = c("a", "a", "b"), taxon = c("x", "y", "x"), quant = c(3, 0, 0))
  got <- as.data.frame(el_frequency_counts(x))

  expect_equal(sorted_frequencies(got), data.frame(sample = "a", count = 3, n_taxa = 1))
})

test_that("non-integer, negative and missing counts are errors that name the column", {
  x <- data.frame(sample = c("a", "a"), taxon = c("x", "y"), n = c(1, 2))
  y <- x
  y$n[1] <- 1.5
  expect_error(el_frequency_counts(y, quant = "n"), "`n`.*integer")
  y$n[1] <- -1
  expect_error(el_frequency_counts(y, quant = "n"), "`n` has negative")
  y$n[1] <- NA
  expect_error(el_frequency_counts(y, quant = "n"), "`n` has missing")
})

test_that("column names are taken from the arguments and the result names its columns", {
  x <- data.frame(site = c("a", "a"), species = c("x", "y"), n = c(1, 1))
  got <- as.data.frame(el_frequency_counts(x, sample = "site", quant = "n"))

  expect_named(got, c("site", "count", "n_taxa"))
  expect_equal(got$n_taxa, 2)
})

test_that("the result is a lazy duckplyr table and the input is converted with prudence 'stingy'", {
  seen <- NULL
  original <- duckplyr::as_duckdb_tibble
  testthat::local_mocked_bindings(
    as_duckdb_tibble = function(x, ..., prudence = "lavish") {
      seen <<- c(seen, prudence)
      original(x, ..., prudence = prudence)
    },
    .package = "duckplyr"
  )

  expect_s3_class(el_frequency_counts(long_table(frequency_fixtures$BCI)), "duckplyr_df")
  expect_true(length(seen) > 0)
  expect_true(all(seen == "stingy"))
})
