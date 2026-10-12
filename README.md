# ecolongr

[![R-CMD-check](https://github.com/biodiversitydata-se/ecolongr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/biodiversitydata-se/ecolongr/actions/workflows/R-CMD-check.yaml)

Community ecology analyses with an aim to complement [Vegan](https://cran.r-project.org/web/packages/vegan/index.html) for large datasets based on long tables.
The priorities are memory economy first and efficiency second.

`ecolongr` computes distances, transformations and ordinations from a table with one row per sample and taxon, without the zeros.
The taxa × samples matrix that `vegan` needs is never built.
Results are tested against `vegan` and have `vegan`'s classes.
`metaMDS()`, `adonis2()`, `betadisper()`, `hclust()` and `cmdscale()` give `vegan`'s results on a distance, and `scores()`, `summary()`, `plot()`, `biplot()`, `envfit()` and `procrustes()` work on an ordination.

## Why

Amplicon data has millions of ASVs and thousands of samples, and almost every cell is zero.
`vegan` needs the dense matrix: 500,000 ASVs × 20,000 samples is about 80 GB of doubles.
`ecolongr` works on the non-zero rows only, in the form the data is delivered in.
For a table of that size, `el_vegdist()` took 25 seconds and 2.4 GB with 12 threads, of which 1.6 GB is the result itself.

## Reading the table

`ecolongr` works on lazy tables from [`duckplyr`](https://duckplyr.tidyverse.org/).
A lazy table is a query that DuckDB runs only when a result is needed, so reading, filtering and joining the table never loads it into R.
Parquet files suit this best, because DuckDB reads only the columns and row groups a query needs.
Text files work too, and keep memory low in the same way:

```r
library(duckplyr)

counts <- read_csv_duckdb(
  "counts.tsv.gz",
  options = list(delim = "\t", header = TRUE),
  prudence = "stingy"
)
```

DuckDB parses a text file again for every query.
To pay that cost once, write the table to Parquet with `compute_parquet(counts, "counts.parquet")` and read that file in the later steps.
A data frame in R is accepted as well, but it is already in memory.

## Status

Early development, version 0.1.0.
The functions below are tested against `vegan`, but names and arguments may still change.
Design, scope and open questions are in [`DESIGN.md`](DESIGN.md).

## Installation

```r
# install.packages("remotes")
remotes::install_github("biodiversitydata-se/ecolongr", build_vignettes = TRUE, dependencies = TRUE)
```

`build_vignettes = TRUE` builds the vignette, which needs `pandoc`; leave it out without `pandoc`, and read [`vignettes/seasonal-ordination.Rmd`](vignettes/seasonal-ordination.Rmd) on GitHub instead.

The package contains C++ code (Rcpp, with OpenMP if your compiler has it), so it needs a C++ compiler.
It uses [`duckplyr`](https://duckplyr.tidyverse.org/) for the lazy table operations.

## Example

The input is a table with a sample, a taxon and a quantity column, with zeros left out.
[nf-core/ampliseq](https://nf-co.re/ampliseq) 3.0.0 and later write such a table to `summary_tables/ampliseq.counts.parquet`.

```r
library(duckplyr)
library(ecolongr)

counts <- read_parquet_duckdb("summary_tables/ampliseq.counts.parquet", prudence = "stingy")
taxonomy <- read_parquet_duckdb(
  "summary_tables/ampliseq.taxonomy.DADA2.consolidated_score.parquet",
  prudence = "stingy"
)

# Select ASVs, rename to the default columns, make relative abundances and compute
# Bray-Curtis distances.
# Nothing is read into R until el_vegdist().
# `samples` sets the order of d to that of the metadata and includes samples without counts.
d <- counts |>
  semi_join(taxonomy, by = "asv_id") |>
  transmute(sample, taxon = asv_id, quant = count) |>
  el_decostand("total") |>
  el_vegdist("bray", samples = metadata$sample)

# d is a dist object: use vegan on it.
# metaMDS() on a dist has no species scores; giving it the matrix would, but needs the dense table.
ordination <- vegan::metaMDS(d)
vegan::adonis2(d ~ treatment, data = metadata)
vegan::envfit(ordination, metadata[, c("temperature", "salinity")], na.rm = TRUE)

# A PCA of the clr of the counts, first five axes.
# `clr` needs a prevalence filter in practice: see "Filter before `clr`" below.
pca <- counts |>
  transmute(sample, taxon = asv_id, quant = count) |>
  el_rda(k = 5, transform = "clr", pseudocount = "inverse_depth", samples = metadata$sample)

summary(pca)
plot(pca)

# Axis labels with the percentage of the total variance: divide by tot.chi, not by sum(eigenvals).
percent <- 100 * vegan::eigenvals(pca)[1:2] / pca$tot.chi
plot(
  vegan::scores(pca, display = "sites", choices = 1:2),
  xlab = sprintf("PC1 (%.1f%%)", percent[1]), ylab = sprintf("PC2 (%.1f%%)", percent[2])
)
vegan::scores(pca, display = "sites", choices = 1:2)
vegan::envfit(pca, metadata[, c("temperature", "salinity")], choices = 1:2, na.rm = TRUE)
```

`ecolongr` does not attach `duckplyr`.
Attach it yourself, as above, to use its verbs on the table.

A longer walk-through, from the published LMO time series to seasonal ordinations, is in the vignette `vignette("seasonal-ordination", package = "ecolongr")`.

### Filter before `clr`

`clr` needs a prevalence filter on amplicon data (see below).
With `duckplyr` verbs the filters stay lazy, so the table is still not read into R:

```r
n_samples <- counts |> summarise(n = n_distinct(sample)) |> pull(n)

# ASVs present in at least 20% of the samples.
common <- counts |>
  summarise(prevalence = n_distinct(sample), .by = asv_id) |>
  filter(prevalence >= 0.2 * n_samples) |>
  select(asv_id)

# Samples with at least 5,000 reads.
deep <- counts |>
  summarise(depth = sum(count, na.rm = TRUE), .by = sample) |>
  filter(depth >= 5000) |>
  select(sample)

# The detection limit of a sample is set by its whole library, so take the depth before filtering.
depth <- counts |>
  summarise(total = sum(count, na.rm = TRUE), .by = sample) |>
  as.data.frame()
pseudocount <- setNames(1 / depth$total, depth$sample)

d <- counts |>
  semi_join(common, by = "asv_id") |>
  semi_join(deep, by = "sample") |>
  transmute(sample, taxon = asv_id, quant = count) |>
  el_vegdist("aitchison", pseudocount = pseudocount)
```

Two details:

- `pseudocount = "inverse_depth"` uses the totals of the table it is given.
  After the filters above those totals are smaller than the libraries, by different amounts in different samples, so the pseudocounts depend on the filter.
  Pass the depth from before filtering by name, as above, unless the table you give is the whole library.
- `duckplyr` translates few R functions.
  `filter(.by = )` is not available (use `summarise(.by = )` and a join as above), `sum()` needs `na.rm = TRUE`, and a scalar with names, such as the result of `quantile()`, must go through `unname()`.

### Data from `asvoccur`

[`asvoccur`](https://github.com/biodiversitydata-se/asvoccur) reads datasets downloaded from the [Swedish ASV portal](https://asv-portal.biodiversitydata.se/).
Download the dataset's ZIP archive from the portal in a browser and put it in a folder, then:

```r
loaded <- asvoccur::load_data("path/to/datasets")
merged <- asvoccur::merge_data(loaded)

# merged$counts is a sparse matrix, taxa in rows and events in columns.
triplets <- Matrix::summary(merged$counts)
counts <- data.frame(
  sample = colnames(merged$counts)[triplets$j],
  taxon = rownames(merged$counts)[triplets$i],
  quant = triplets$x
)
d <- el_vegdist(counts, "aitchison", pseudocount = "inverse_depth")
```

`merged$events` has the sample metadata, with `eventID` matching the sample names.
`asvoccur` does not yet have a long-table form of its own; [asvoccur#1](https://github.com/biodiversitydata-se/asvoccur/issues/1) proposes one.

## Functions

| function | what it does | `vegan` equivalent |
|---|---|---|
| `el_vegdist()` | `bray`, `jaccard`, `euclidean`, `manhattan`, `canberra`, `aitchison` and `robust.aitchison` dissimilarities, optionally `binary`; returns a `dist` | `vegdist()` |
| `el_decostand()` | `total`, `hellinger`, `log`, `pa` and `rclr`, as a lazy long table | `decostand()` |
| `el_rda()` | unconstrained PCA of the table, of its clr or of its completed rclr, the first `k` axes; returns an `rda` object | `rda()` |

`metaMDS()` and `adonis2()` only need a `dist`, so use `vegan`'s own on the result of `el_vegdist()`.

`clr` is not a method of `el_decostand()`, because it makes every cell non-zero, which a table without zeros cannot hold.
It is computed from the raw counts inside the two functions that use it.
`el_vegdist(method = "aitchison")` is the Euclidean distance between clr-transformed samples, which is the Aitchison distance.
`el_rda(transform = "clr")` is a PCA of the clr, and with all axes its sample distances are the same Aitchison distances.
The `pseudocount` is in the units of the counts: a positive number, `"inverse_depth"` for `1 / (the sample's total count)`, or a numeric vector named by sample.
`"inverse_depth"` uses the totals of the table you pass, so on a filtered table pass the depths from before filtering as a named vector (see "Filter before `clr`").

Conventions:

- The default column names are `sample`, `taxon` and `quant`.
  Each function has arguments with the same names to point at other columns.
- Rows with `quant == 0` are read as absent cells, but their sample still counts: a sample with only such rows is an all-zero sample.
- A sample with no rows at all cannot be seen in the table.
  Pass `samples = ` with all sample IDs to include them as all-zero samples, as `vegan` does with a zero row.
  Without it, the samples that occur in the table are used, in sorted order.
  Pass `samples = metadata$sample` to `el_vegdist()` and `el_rda()` before using `envfit()`, `adonis2()` or `procrustes()` with the metadata, which pair rows by position.
- Bray–Curtis and Jaccard need non-negative quantities.
  Euclidean distance and `el_rda()` accept negative ones, for example after `rclr`.

## How results relate to `vegan`

Distances and standardisations equal `vegan`'s to about 1e-10, and PCA scores to about 1e-8 (1e-6 with the iterative solver on large tables), apart from the cases below.
The tests compare every function and option with `vegan` on `varespec` and a synthetic community, including all-zero samples and explicit zero rows.

- **`el_decostand("rclr")` leaves zeros absent.**
  `vegan`'s `decostand(x, "rclr")` fills the zeros with values from a matrix completion, so every cell becomes non-zero.
  Like `clr`, that cannot be returned as a table without zeros.
  The completion is also one global fit over all observed cells, which a lazy table query cannot express.
  `el_decostand("rclr")` therefore does only the part that fits: the log-ratio of the non-zero cells.
  The completion also changes the non-zero cells; `el_decostand("rclr")` matches `decostand(x, "rclr", impute = FALSE)`.
  A later step that reads absent cells as 0 puts those taxa at the sample's geometric mean, so Euclidean distance on it is not `vegan`'s `robust.aitchison` when the table has zeros.
  For that, use `el_vegdist(method = "robust.aitchison")` or `el_rda(transform = "rclr")` on the raw counts: they run the same matrix completion as `vegan`, of rank 3 unless you set `ropt`, on the observed cells only and agree with `vegan` to about 1e-6.
  They need every sample to have at least one non-zero cell, and `el_rda(transform = "rclr")` returns at most `ropt` axes (three by default, as in `vegan`).
- **`clr` needs a prevalence filter in practice.**
  On a 414-sample, 60,633-ASV time series with 99% zeros, PC1–PC2 of a `clr` PCA showed the seasonal cycle poorly (seasonal R² 0.28), but as well as Hellinger once ASVs present in fewer than 20% of the samples were removed (0.64).
  Filter with `duckplyr` verbs before the call, as in "Filter before `clr`" above: `semi_join()` against the ASVs with `n_distinct(sample) >= 0.2 * n_distinct(all samples)`.
- **`metaMDS()` on the `dist` has no species scores.**
  Given a community matrix, `metaMDS()` also computes species scores as weighted averages and applies its own automatic transformations (`autotransform`).
  That needs the dense matrix, which is what this package avoids.
  On a `dist` you get the sample scores only: `ordination$species` is `NA` and `scores(ordination, display = "species")` is empty.
  Transform with `el_decostand()` before `el_vegdist()` instead.
- **`log` follows `vegan` without the warning:** if any value is not a whole number, everything is first divided by the smallest positive value in the table.
- **`el_rda()` returns the first `k` axes only.**
  `scores()`, `eigenvals()`, `plot()` and `screeplot()` work.
  `predict()` and `fitted()` need the full matrix and give an error.
  `anova()` gives `vegan`'s table for an unconstrained ordination, which has no constrained component.
  Axis signs are arbitrary, as in `vegan`.
- **`k` sets the number of axes (default 5).**
  A `k` above the rank, `min(number of samples − 1, number of taxa)`, is silently reduced to it, so a very large `k` returns all axes, as a full PCA.
  A `k` close to the rank makes `el_rda()` use the exact solver, which builds the dense samples × taxa matrix.
  Keep `k` small on large tables, where the iterative solver is used.
  With `transform = "rclr"` there are at most `ropt` axes, the rank of the completion (default 3, as `vegan`'s `ropt`).
- **The variance the `k` axes leave out is known, and `el_rda()` reports it in a message.**
  `pca$tot.chi` is the total variance of the whole table, whatever `k` is, so `1 - sum(vegan::eigenvals(pca)) / pca$tot.chi` is the unexplained share.
  Only that sum is known: the eigenvalues of the axes beyond `k` are not.
  Silence the message with `suppressMessages()` or `options(ecolongr.quiet = TRUE)`.
- **Divide an eigenvalue by `tot.chi`, never by the sum of the eigenvalues.**
  `vegan::eigenvals(pca) / sum(vegan::eigenvals(pca))` sums to one over the `k` axes whatever `k` is, and so hides the variance that the other axes hold.
  `vegan::eigenvals(pca) / pca$tot.chi` and `summary(pca)` give the share of the whole.
- **Some plotting packages divide by the sum of the eigenvalues and so mislabel the axes.**
  `ggordiplots` (`ord_labels()`, used by its ordination plot functions) and `microViz` (`ord_plot()`) put a percentage in the axis label that is the eigenvalue divided by the sum of the kept eigenvalues, so the `k` labels sum to 100%.
  Labels from `eigenvals(pca) / pca$tot.chi`, as above, are right; replace the axis titles of such a plot with them, for example with `ggplot2::xlab()`.
  `ggvegan::autoplot()` does not print percentages and works on the result.
  For `transform = "rclr"` the share is 0, because the completed matrix has rank `ropt`, and the axes hold all of it.
- **`hellinger` writes a temporary Parquet file**, because `duckplyr` cannot compute a square root.
  The file is removed when the R session ends.

## Memory and speed

`el_vegdist()` runs the lazy query in DuckDB, writes it to a temporary Parquet file sorted by taxon, and reads it back in batches of whole taxa.
A compiled kernel adds up the pairs into one vector in the order of a `dist` object.
The memory has three parts:

- the result, one number per pair of samples, which grows with the square of the number of samples;
- one batch of rows;
- the memory DuckDB uses for the sort, which `memory_limit` caps (250 MB by default).

Measured on a 24-thread desktop under load, on synthetic tables where a few taxa occur in nearly every sample:

| table | `el_vegdist()` | `el_rda()` |
|---|---|---|
| 5,000 samples, 11.1M rows | 7 s with 12 threads (30 s with one), 1.0 GB | 35 s, 1.4 GB |
| 20,000 samples, 19.4M rows | 25 s with 12 threads (148 s with one), 2.4 GB | 63 s, 2.2 GB |

Settings:

- `threads` (option `ecolongr.threads`, default half the cores) sets the threads of the distance kernel.
  The result is identical for any number.
- `memory_limit` (option `ecolongr.memory_limit`, default `"250MB"`) limits the memory DuckDB uses for sorting, and it spills to disk beyond that.
- `batch_rows` (default 1,000,000) is roughly how many table rows are in R at a time.

`ecolongr` does not change the settings of your `duckplyr` session.
The lazy steps you write before calling it, such as `el_decostand()`, run there with DuckDB's defaults (all cores, 80% of the memory).
On a small machine, limit them:

```r
duckplyr::db_exec("SET memory_limit = '1GB'")
duckplyr::db_exec("SET threads = 4")
```

On the 5,000-sample table this took the peak of a `rclr` plus Euclidean pipeline from 2.4 GB to 1.1 GB at a cost of 3 seconds.

Make the first table `prudence = "stingy"`, as in the example, to have `duckplyr` raise an error where it would otherwise silently fall back to `dplyr` and copy the whole table into R.

## Development

- [`DESIGN.md`](DESIGN.md) has the scope, decisions, measurements and open questions.
- [`AGENTS.md`](AGENTS.md) has the working rules: tests first, `vegan` as the reference.
- `bench/` has scripts for the large synthetic tables.
  Benchmark with an installed package, not `devtools::load_all()`, which compiles without optimisation.
- `devtools::test()` runs the tests.

## License

MIT, see `LICENSE.md`.
