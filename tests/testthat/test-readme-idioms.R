# The prevalence and depth filters shown in the README must stay lazy on a stingy
# duckplyr table, and must select what a plain R computation selects.
library(dplyr, warn.conflicts = FALSE)

lmo <- lmo_subset()
counts <- duckplyr::as_duckdb_tibble(
  data.frame(taxon = lmo$long$taxon, sample = lmo$long$sample, quant = lmo$long$quant),
  prudence = "stingy"
)

test_that("a prevalence filter written with duckplyr verbs stays lazy and matches base R", {
  n_samples <- counts |> summarise(n = n_distinct(sample)) |> pull(n)
  common <- counts |>
    summarise(prevalence = n_distinct(sample), .by = taxon) |>
    filter(prevalence >= 0.2 * n_samples) |>
    select(taxon)
  kept <- counts |> semi_join(common, by = "taxon")

  expect_s3_class(kept, "duckplyr_df")
  expected <- colnames(lmo$matrix)[colSums(lmo$matrix > 0) >= 0.2 * nrow(lmo$matrix)]
  expect_setequal(as.data.frame(kept)$taxon, expected)
  expect_equal(nrow(as.data.frame(kept)), sum(lmo$matrix[, expected] > 0))
})

test_that("a depth filter written with duckplyr verbs stays lazy and matches base R", {
  threshold <- unname(stats::quantile(rowSums(lmo$matrix), 0.2))
  deep <- counts |>
    summarise(depth = sum(quant, na.rm = TRUE), .by = sample) |>
    filter(depth >= threshold) |>
    select(sample)
  kept <- counts |> semi_join(deep, by = "sample")

  expect_s3_class(kept, "duckplyr_df")
  expect_setequal(unique(as.data.frame(kept)$sample), rownames(lmo$matrix)[rowSums(lmo$matrix) >= threshold])
})

test_that("the per-sample depth can be taken before filtering and passed as a named pseudocount", {
  depth <- counts |> summarise(total = sum(quant, na.rm = TRUE), .by = sample) |> as.data.frame()
  pseudocount <- stats::setNames(1 / depth$total, depth$sample)
  common <- counts |>
    summarise(prevalence = n_distinct(sample), .by = taxon) |>
    filter(prevalence >= 6) |>
    select(taxon)
  d <- counts |>
    semi_join(common, by = "taxon") |>
    select(sample, taxon, quant) |>
    el_vegdist("aitchison", pseudocount = pseudocount)

  m <- lmo$matrix[, colSums(lmo$matrix > 0) >= 6]
  expect_equal(as_full_matrix(d), as_full_matrix(stats::dist(clr_reference(m, pseudocount))), tolerance = 1e-10)
})

test_that("inverse_depth uses the totals of the table it is given, which is the filtered one", {
  common <- counts |>
    summarise(prevalence = n_distinct(sample), .by = taxon) |>
    filter(prevalence >= 6) |>
    select(taxon)
  filtered <- counts |> semi_join(common, by = "taxon") |> select(sample, taxon, quant)
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

test_that("axis labels with the percentage of the whole variance equal summary()'s proportions, for k below the rank", {
  m <- present_taxa(lmo$matrix)
  pca <- el_rda(lmo$long, k = 3, transform = "clr", pseudocount = 1)
  percent <- 100 * vegan::eigenvals(pca)[1:2] / pca$tot.chi
  ref <- 100 * summary(vegan::rda(clr_reference(lmo$matrix[, colSums(lmo$matrix) > 0], 1)))$cont$importance["Proportion Explained", 1:2]

  expect_equal(as.numeric(percent), unname(ref), tolerance = 1e-6)
  expect_lt(sum(percent), 99)
  expect_match(sprintf("PC1 (%.1f%%)", percent[1]), "^PC1 \\(\\d+\\.\\d%\\)$")
})

test_that("the alpha diversity idioms in the README line up by sample name and match vegan", {
  richness <- el_specnumber(counts)
  shannon <- el_diversity(counts, "shannon")
  diversity <- tibble(sample = names(shannon), shannon = shannon, richness = richness[sample])

  expect_equal(diversity$richness, vegan::specnumber(lmo$matrix)[diversity$sample], ignore_attr = TRUE)
  expect_equal(diversity$shannon, vegan::diversity(lmo$matrix, "shannon")[diversity$sample], ignore_attr = TRUE)

  depth <- counts |>
    summarise(reads = sum(quant, na.rm = TRUE), .by = sample) |>
    summarise(depth = min(reads)) |>
    pull(depth)
  expect_equal(depth, min(rowSums(lmo$matrix)))

  expected_richness <- el_rarefy(counts, depth)
  expect_equal(expected_richness, vegan::rarefy(lmo$matrix, depth)[names(expected_richness)], ignore_attr = TRUE, tolerance = 1e-10)
  expect_true(all(expected_richness <= richness[names(expected_richness)] + 1e-9))

  set.seed(1)
  subsample <- as.data.frame(el_rrarefy(counts))
  expect_true(all(tapply(subsample$quant, subsample$sample, sum) == depth))
})

test_that("the frequency table splits into the two-column tables that breakaway reads", {
  frequencies <- counts |> el_frequency_counts() |> as.data.frame()
  per_sample <- frequencies |>
    split(frequencies$sample) |>
    purrr::map(\(d) select(d, count, n_taxa))

  expect_equal(length(per_sample), nrow(lmo$matrix))
  expect_true(all(purrr::map_lgl(per_sample, \(d) identical(names(d), c("count", "n_taxa")))))
  expect_equal(sum(per_sample[[1]]$count * per_sample[[1]]$n_taxa), sum(lmo$matrix[names(per_sample)[1], ]))
})
