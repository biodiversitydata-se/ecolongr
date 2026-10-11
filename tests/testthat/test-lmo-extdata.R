lmo_files <- function() {
  list(
    counts = system.file("extdata", "lmo_free_living.tsv.gz", package = "ecolongr"),
    samples = system.file("extdata", "lmo_free_living_samples.tsv", package = "ecolongr")
  )
}

test_that("the bundled LMO counts are a long table without zeros, with one row per sample and taxon", {
  files <- lmo_files()
  expect_true(nzchar(files$counts))
  counts <- utils::read.delim(gzfile(files$counts), stringsAsFactors = FALSE)

  expect_named(counts, c("sample", "taxon", "quant"))
  expect_true(all(counts$quant >= 1))
  expect_equal(counts$quant, round(counts$quant))
  expect_false(anyDuplicated(counts[c("sample", "taxon")]) > 0)
  expect_false(anyNA(counts))
})

test_that("every bundled taxon occurs in at least 5% of the bundled samples", {
  files <- lmo_files()
  counts <- utils::read.delim(gzfile(files$counts), stringsAsFactors = FALSE)
  samples <- utils::read.delim(files$samples, stringsAsFactors = FALSE)

  prevalence <- tapply(counts$sample, counts$taxon, function(s) length(unique(s)))
  expect_gte(min(prevalence), 0.05 * nrow(samples))
  expect_gt(max(prevalence), 0.9 * nrow(samples))
})

test_that("the sample table covers exactly the samples of the counts and has dates and a day of year", {
  files <- lmo_files()
  counts <- utils::read.delim(gzfile(files$counts), stringsAsFactors = FALSE)
  samples <- utils::read.delim(files$samples, stringsAsFactors = FALSE)

  expect_named(samples, c("sample", "date", "depth", "temperature", "salinity"))
  expect_setequal(samples$sample, counts$sample)
  expect_false(anyDuplicated(samples$sample) > 0)
  dates <- as.Date(samples$date)
  expect_false(anyNA(dates))
  expect_gte(min(dates), as.Date("2011-01-01"))
  expect_lte(max(dates), as.Date("2023-12-31"))
  expect_gt(length(unique(format(dates, "%Y"))), 10)
})

test_that("the samples span the seasons evenly enough for a seasonal analysis", {
  samples <- utils::read.delim(lmo_files()$samples, stringsAsFactors = FALSE)
  months <- table(factor(format(as.Date(samples$date), "%m"), levels = sprintf("%02d", 1:12)))

  expect_true(all(months >= 5))
})

test_that("temperature is complete and salinity mostly so, both in plausible units for the Baltic Sea", {
  samples <- utils::read.delim(lmo_files()$samples, stringsAsFactors = FALSE)

  expect_type(samples$temperature, "double")
  expect_false(anyNA(samples$temperature))
  expect_true(all(samples$temperature > -2 & samples$temperature < 30))
  expect_type(samples$salinity, "double")
  expect_true(all(samples$salinity[!is.na(samples$salinity)] > 3 & samples$salinity[!is.na(samples$salinity)] < 12))
  expect_gt(sum(is.na(samples$salinity)), 0)
  expect_lt(mean(is.na(samples$salinity)), 0.1)
})

test_that("temperature follows the seasons: warmest in summer, coldest in winter", {
  samples <- utils::read.delim(lmo_files()$samples, stringsAsFactors = FALSE)
  month <- as.integer(format(as.Date(samples$date), "%m"))
  by_month <- tapply(samples$temperature, month, mean)

  expect_true(names(which.max(by_month)) %in% c("7", "8"))
  expect_true(names(which.min(by_month)) %in% c("1", "2", "3"))
})
