# The prevalence and depth filters shown in the README must stay lazy on a stingy
# duckplyr table, and must select what a plain R computation selects.
library(dplyr, warn.conflicts = FALSE)

lmo <- lmo_subset()
counts <- duckplyr::as_duckdb_tibble(
  data.frame(asv_id = lmo$long$taxon, sample = lmo$long$sample, count = lmo$long$quant),
  prudence = "stingy"
)

test_that("a prevalence filter written with duckplyr verbs stays lazy and matches base R", {
  n_samples <- counts |> summarise(n = n_distinct(sample)) |> pull(n)
  common <- counts |>
    summarise(prevalence = n_distinct(sample), .by = asv_id) |>
    filter(prevalence >= 0.2 * n_samples) |>
    select(asv_id)
  kept <- counts |> semi_join(common, by = "asv_id")

  expect_s3_class(kept, "duckplyr_df")
  expected <- colnames(lmo$matrix)[colSums(lmo$matrix > 0) >= 0.2 * nrow(lmo$matrix)]
  expect_setequal(as.data.frame(kept)$asv_id, expected)
  expect_equal(nrow(as.data.frame(kept)), sum(lmo$matrix[, expected] > 0))
})

test_that("a depth filter written with duckplyr verbs stays lazy and matches base R", {
  threshold <- unname(stats::quantile(rowSums(lmo$matrix), 0.2))
  deep <- counts |>
    summarise(depth = sum(count, na.rm = TRUE), .by = sample) |>
    filter(depth >= threshold) |>
    select(sample)
  kept <- counts |> semi_join(deep, by = "sample")

  expect_s3_class(kept, "duckplyr_df")
  expect_setequal(unique(as.data.frame(kept)$sample), rownames(lmo$matrix)[rowSums(lmo$matrix) >= threshold])
})

test_that("the per-sample depth can be taken before filtering and passed as a named pseudocount", {
  depth <- counts |> summarise(total = sum(count, na.rm = TRUE), .by = sample) |> as.data.frame()
  pseudocount <- stats::setNames(1 / depth$total, depth$sample)
  common <- counts |>
    summarise(prevalence = n_distinct(sample), .by = asv_id) |>
    filter(prevalence >= 6) |>
    select(asv_id)
  d <- counts |>
    semi_join(common, by = "asv_id") |>
    transmute(sample, taxon = asv_id, quant = count) |>
    el_vegdist("aitchison", pseudocount = pseudocount)

  m <- lmo$matrix[, colSums(lmo$matrix > 0) >= 6]
  expect_equal(as_full_matrix(d), as_full_matrix(stats::dist(clr_reference(m, pseudocount))), tolerance = 1e-10)
})

test_that("inverse_depth uses the totals of the table it is given, which is the filtered one", {
  common <- counts |>
    summarise(prevalence = n_distinct(sample), .by = asv_id) |>
    filter(prevalence >= 6) |>
    select(asv_id)
  filtered <- counts |> semi_join(common, by = "asv_id") |> transmute(sample, taxon = asv_id, quant = count)
  m <- lmo$matrix[, colSums(lmo$matrix > 0) >= 6]

  d <- el_vegdist(filtered, "aitchison", pseudocount = "inverse_depth")

  expect_equal(as_full_matrix(d), as_full_matrix(stats::dist(clr_reference(m, 1 / rowSums(m)))), tolerance = 1e-10)
})

test_that("a taxa-by-events sparse matrix, as asvoccur returns, converts to a long table", {
  m <- lmo$matrix
  counts_matrix <- Matrix::Matrix(t(m), sparse = TRUE)
  triplets <- Matrix::summary(counts_matrix)
  long <- data.frame(sample = colnames(counts_matrix)[triplets$j], taxon = rownames(counts_matrix)[triplets$i], quant = triplets$x)

  expect_equal(nrow(long), sum(m > 0))
  expect_equal(as_full_matrix(el_vegdist(long, "bray")), as_full_matrix(vegan::vegdist(m, "bray")), tolerance = 1e-10)
})
