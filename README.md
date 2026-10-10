# ecolongr

Community ecology analyses on long tables.

`ecolongr` computes distances, transformations and ordinations from a table with one row per sample and taxon, without the zeros.
The taxa × samples matrix that `vegan` needs is never built.
Results are tested against `vegan` and have `vegan`'s classes, so `vegan::metaMDS()`, `vegan::adonis2()`, `scores()` and `plot()` work on them.

## Why

Amplicon data has millions of ASVs and thousands of samples, and almost every cell is zero.
`vegan` needs the dense matrix: 500,000 ASVs × 20,000 samples is about 80 GB of doubles.
`ecolongr` works on the non-zero rows only, in the form the data is delivered in.
For a table of that size, `el_vegdist()` took 25 seconds and 2.4 GB with 12 threads, of which 1.6 GB is the result itself.

## Status

Early development, version 0.0.0.9000.
The functions below are tested against `vegan`, but names and arguments may still change.
Design, scope and open questions are in [`DESIGN.md`](DESIGN.md).

## Installation

```r
# install.packages("remotes")
remotes::install_github("biodiversitydata-se/ecolongr")
```

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
  "summary_tables/ampliseq.taxonomy.DADA2.consolidated_most-specific.parquet",
  prudence = "stingy"
)

# Select ASVs, rename to the default columns, standardise and compute distances.
# Nothing is read into R until el_vegdist().
# `samples` sets the order of d to that of the metadata and includes samples without counts.
d <- counts |>
  semi_join(taxonomy, by = "asv_id") |>
  transmute(sample, taxon = asv_id, quant = count) |>
  el_decostand("hellinger") |>
  el_vegdist("euclidean", samples = metadata$sample)

# d is a dist object: use vegan on it.
ordination <- vegan::metaMDS(d)
vegan::adonis2(d ~ treatment, data = metadata)

# A PCA of the same table, first five axes.
pca <- counts |>
  transmute(sample, taxon = asv_id, quant = count) |>
  el_decostand("rclr") |>
  el_rda(k = 5)

vegan::eigenvals(pca) / pca$tot.chi
plot(pca)
vegan::scores(pca, display = "sites", choices = 1:2)
```

`ecolongr` does not attach `duckplyr`.
Attach it yourself, as above, to use its verbs on the table.

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

d <- counts |>
  semi_join(common, by = "asv_id") |>
  semi_join(deep, by = "sample") |>
  transmute(sample, taxon = asv_id, quant = count) |>
  el_vegdist("aitchison", pseudocount = "inverse_depth")
```

Two details:

- `pseudocount = "inverse_depth"` uses the totals of the table it is given, which here is the filtered one.
  To use each sample's depth from before filtering, compute it first and pass it by name:
  `depth <- counts |> summarise(total = sum(count, na.rm = TRUE), .by = sample) |> as.data.frame()`, then
  `pseudocount = setNames(1 / depth$total, depth$sample)`.
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
| `el_vegdist()` | `bray`, `jaccard`, `euclidean`, `manhattan`, `canberra` and `aitchison` dissimilarities, optionally `binary`; returns a `dist` | `vegdist()` |
| `el_decostand()` | `total`, `hellinger`, `log`, `pa` and `rclr`, as a lazy long table | `decostand()` |
| `el_rda()` | unconstrained PCA of the table or of its clr, the first `k` axes; returns an `rda` object | `rda()` |

`metaMDS()` and `adonis2()` only need a `dist`, so use `vegan`'s own on the result of `el_vegdist()`.

`clr` is not a method of `el_decostand()`, because it makes every cell non-zero.
Use `el_vegdist(method = "aitchison")` or `el_rda(transform = "clr")` on the raw counts instead.
The `pseudocount` is in the units of the counts: a positive number, `"inverse_depth"` for `1 / (the sample's total count)`, or a numeric vector named by sample.

Conventions:

- The default column names are `sample`, `taxon` and `quant`.
  Each function has arguments with the same names to point at other columns.
- Rows with `quant == 0` are ignored.
- A sample with no non-zero rows cannot be seen in the table.
  Pass `samples = ` with all sample IDs to include them as all-zero samples, as `vegan` does with a zero row.
  Without it, the samples that occur in the table are used.
- Bray–Curtis and Jaccard need non-negative quantities.
  Euclidean distance and `el_rda()` accept negative ones, for example after `rclr`.

## How results relate to `vegan`

Distances and standardisations equal `vegan`'s to about 1e-10, and PCA scores to about 1e-8 (1e-6 with the iterative solver on large tables), apart from the cases below.
The tests compare every function and option with `vegan` on `varespec` and a synthetic community, including all-zero samples and explicit zero rows.

- **`rclr` leaves zeros absent.**
  `vegan`'s `decostand(x, "rclr")` fills them by matrix completion, which is a global step that also changes the non-zero cells.
  `el_decostand("rclr")` matches `decostand(x, "rclr", impute = FALSE)`.
  A later step that reads absent cells as 0 puts those taxa at the sample's geometric mean.
  Euclidean distance on that is therefore not `vegan`'s `robust.aitchison` when the table has zeros, but it is identical when it has none.
- **`clr` needs a prevalence filter in practice.**
  On a 414-sample, 60,633-ASV time series with 99% zeros, PC1–PC2 of a `clr` PCA showed the seasonal cycle poorly (seasonal R² 0.28), but as well as Hellinger once ASVs present in fewer than 20% of the samples were removed (0.64).
  Filter with `duckplyr` verbs before the call, for example `semi_join()` against the ASVs with `n_distinct(sample) >= 0.2 * n_distinct(all samples)`.
- **`log` follows `vegan` without the warning:** if any value is not a whole number, everything is first divided by the smallest positive value in the table.
- **`el_rda()` returns the first `k` axes only.**
  `scores()`, `eigenvals()`, `plot()` and `screeplot()` work.
  `predict()`, `anova()` and `fitted()` need the full matrix and do not.
  Axis signs are arbitrary, as in `vegan`.
- **`hellinger` writes a temporary Parquet file**, because `duckplyr` cannot compute a square root.
  The file is removed when the R session ends.

## Memory and speed

`el_vegdist()` runs the lazy query in DuckDB, writes it to a temporary Parquet file sorted by taxon, and reads it back in batches of whole taxa.
A compiled kernel adds up the pairs into one vector in the order of a `dist` object.
The memory is the result (one number per pair of samples, so it grows with the square of the number of samples), one batch, and DuckDB's sort buffer.

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
