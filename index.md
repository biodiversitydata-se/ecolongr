# ecolongr

[![R-CMD-check](https://github.com/biodiversitydata-se/ecolongr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/biodiversitydata-se/ecolongr/actions/workflows/R-CMD-check.yaml)

Community ecology analyses with an aim to complement
[Vegan](https://cran.r-project.org/web/packages/vegan/index.html) for
large datasets based on long tables. The priorities are memory economy
first and efficiency second.

`ecolongr` computes diversity, distances, transformations and
ordinations from a table with one row per sample and taxon, without the
zeros. The taxa × samples matrix that `vegan` needs is never built.
Results are tested against `vegan` and have `vegan`’s classes.
`metaMDS()`, `adonis2()`, `betadisper()`,
[`hclust()`](https://rdrr.io/r/stats/hclust.html) and
[`cmdscale()`](https://rdrr.io/r/stats/cmdscale.html) give `vegan`’s
results on a distance, and `scores()`,
[`summary()`](https://rdrr.io/r/base/summary.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html),
[`biplot()`](https://rdrr.io/r/stats/biplot.html), `envfit()` and
`procrustes()` work on an ordination.

The functions are documented one by one on the [package
website](https://biodiversitydata-se.github.io/ecolongr/reference/index.html),
with examples that run on the table bundled in the package.

## Why

Amplicon data has millions of ASVs and thousands of samples, and almost
every cell is zero. `vegan` needs the dense matrix: 500,000 ASVs ×
20,000 samples is about 80 GB of doubles. `ecolongr` works on the
non-zero rows only, in the form the data is delivered in. For a table of
that size,
[`el_vegdist()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_vegdist.md)
took 25 seconds and 2.4 GB with 12 threads, of which 1.6 GB is the
result itself.

Bray-Curtis distances on synthetic tables, one run each; the memory is
the peak of the whole R process:

| table | `ecolongr` | `vegan` |
|----|----|----|
| 500 samples, 41,090 taxa | 0.4 s, 270 MB | 47 s, 530 MB |
| 2,000 samples, 158,708 taxa | 3.1 s, 390 MB | more than 33 min, at least 4.1 GB |
| 3,000 samples, 235,740 taxa | 5.7 s, 480 MB | not run |

`ecolongr` ran on one thread, as `vegan` does. The slow `vegan` run on
the 2,000-sample table was stopped by hand, so its memory is what it
held at that point. The 3,000-sample table was run only with `ecolongr`:
scaling the 500-sample time with the square of the samples and the
number of taxa, `vegan` would need about 2.5 hours. The timings of
`vegan` varied by a factor of two between runs on a machine under load.
The Bioconductor package [`mia`](https://bioconductor.org/packages/mia)
calls
[`vegan::vegdist()`](https://vegandevs.github.io/vegan/reference/vegdist.html)
on the dense matrix, so it takes as long as `vegan` does.

## Status

Early development, version 0.1.1. The functions are tested against
`vegan`, but names and arguments may still change. Design, scope and
open questions are in
[`DESIGN.md`](https://biodiversitydata-se.github.io/ecolongr/DESIGN.md).

## Installation

``` r

# install.packages("remotes")
remotes::install_github("biodiversitydata-se/ecolongr", build_vignettes = TRUE, dependencies = TRUE)
```

`build_vignettes = TRUE` builds the vignettes, which needs `pandoc`;
leave it out without `pandoc`, and read them on the [package
website](https://biodiversitydata-se.github.io/ecolongr/articles/index.html)
instead.

The package contains C++ code (Rcpp, with OpenMP if your compiler has
it), so it needs a C++ compiler. It uses
[`duckplyr`](https://duckplyr.tidyverse.org/) for the lazy table
operations.

## Reading a table

The input is a table with a sample, a taxon and a quantity column, with
zeros left out. `ecolongr` works on lazy tables from
[`duckplyr`](https://duckplyr.tidyverse.org/). A lazy table is a query
that DuckDB runs only when a result is needed, so reading, filtering and
joining the table never loads it into R. Parquet files suit this best,
because DuckDB reads only the columns and row groups a query needs.
[nf-core/ampliseq](https://nf-co.re/ampliseq) 3.0.0 and later write such
a table to `summary_tables/ampliseq.counts.parquet`.

``` r

library(duckplyr)
library(ecolongr)

counts <- read_parquet_duckdb("summary_tables/ampliseq.counts.parquet", prudence = "stingy") |>
  transmute(sample, taxon = asv_id, quant = count)
```

`prudence = "stingy"` makes `duckplyr` stop with an error where it would
otherwise silently copy the whole table into R. Text files work too, and
keep memory low in the same way:

``` r

counts <- read_csv_duckdb(
  "counts.tsv.gz",
  options = list(delim = "\t", header = TRUE),
  prudence = "stingy"
)
```

DuckDB parses a text file again for every query. To pay that cost once,
write the table to Parquet with
`compute_parquet(counts, "counts.parquet")` and read that file in the
later steps. A data frame in R is accepted as well, but it is already in
memory. `ecolongr` does not attach `duckplyr`: attach it yourself, as
above, to use its verbs on the table.

The default column names are `sample`, `taxon` and `quant`, and every
function has arguments with the same names to point at other columns.
The package bundles a small table, 347 samples of a 16S time series from
the Baltic Sea, that the examples in the function help pages use:

``` r

counts <- read_csv_duckdb(
  system.file("extdata", "lmo_free_living.tsv.gz", package = "ecolongr"),
  options = list(delim = "\t", header = TRUE), prudence = "stingy"
)
metadata <- read.delim(system.file("extdata", "lmo_free_living_samples.tsv", package = "ecolongr"))
```

A table without zeros cannot show a sample that has no rows at all. Pass
`samples = metadata$sample` to
[`el_vegdist()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_vegdist.md)
and
[`el_rda()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_rda.md)
to include such samples, as `vegan` does with a zero row, and to set the
order of the result to that of the metadata. This matters before
`envfit()`, `adonis2()` or `procrustes()`, which pair rows by position.

### Data from `asvoccur`

[`asvoccur`](https://github.com/biodiversitydata-se/asvoccur) reads
datasets downloaded from the [Swedish ASV
portal](https://asv-portal.biodiversitydata.se/). It returns a sparse
matrix, which `Matrix::summary()` turns into a long table:

``` r

loaded <- asvoccur::load_data("path/to/datasets")
merged <- asvoccur::merge_data(loaded)

triplets <- Matrix::summary(merged$counts)
counts <- data.frame(
  sample = colnames(merged$counts)[triplets$j],
  taxon = rownames(merged$counts)[triplets$i],
  quant = triplets$x
)
```

`merged$events` has the sample metadata, with `eventID` matching the
sample names. `asvoccur` does not yet have a long-table form of its own;
[asvoccur#1](https://github.com/biodiversitydata-se/asvoccur/issues/1)
proposes one.

## Alpha diversity

Richness and diversity indices of each sample
([`el_specnumber()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_diversity.html),
[`el_diversity()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_diversity.html)),
richness at a common depth
([`el_rarefy()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_drarefy.html),
[`el_drarefy()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_drarefy.html),
[`el_rrarefy()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_rrarefy.html))
and the frequency table that richness estimators need
([`el_frequency_counts()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_frequency_counts.html)).
They correspond to
[`vegan::specnumber()`](https://vegandevs.github.io/vegan/reference/diversity.html),
`diversity()`, `rarefy()`, `drarefy()` and `rrarefy()`.

The first two return one number per sample, as a vector named by sample:

``` r

richness <- el_specnumber(counts)
shannon <- el_diversity(counts, "shannon")      # or "simpson", "invsimpson"; `base` sets the logarithm

diversity <- tibble(sample = names(shannon), shannon = shannon, richness = richness[sample]) |>
  left_join(metadata, by = "sample")
```

Richness grows with the number of reads, so compare samples at a common
depth.
[`el_rarefy()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_drarefy.md)
gives the richness each sample is expected to have in a subsample of
`depth` reads, which is exact and needs no random numbers.
[`el_rrarefy()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_rrarefy.md)
draws one random subsample, as
[`vegan::rrarefy()`](https://vegandevs.github.io/vegan/reference/rarefy.html)
does, and returns a table; use it when you need the subsampled table
itself.

``` r

depth <- counts |>
  summarise(reads = sum(quant, na.rm = TRUE), .by = sample) |>
  summarise(depth = min(reads)) |>
  pull(depth)

expected_richness <- el_rarefy(counts, depth)

set.seed(1)
subsample <- el_rrarefy(counts, depth)
```

Estimators of the total richness, including the taxa not seen, such as
[`breakaway`](https://adw96.github.io/breakaway/), fit a model to the
number of taxa seen once, twice and so on.
[`el_frequency_counts()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_frequency_counts.md)
computes that table for all samples in one lazy step, and the fits then
run on small tables, one per sample, in parallel if you like:

``` r

frequencies <- counts |> el_frequency_counts() |> as.data.frame()

# breakaway reads one table per sample, sorted by count.
per_sample <- frequencies |>
  arrange(sample, count) |>
  split(~ sample) |>
  purrr::map(\(d) select(d, count, n_taxa))

estimates <- per_sample |>
  purrr::map(purrr::in_parallel(\(d) breakaway::breakaway(d)))   # after mirai::daemons(4)
```

`breakaway`’s estimates depend strongly on its `cutoff` argument, and
tables without singletons, such as DADA2 output, need
`breakaway::breakaway_nof1()`. `ecolongr` does not reimplement these
estimators.

## Beta diversity and ordination

The steps from a table to a distance matrix and an ordination:

1.  Filter the table with `duckplyr` verbs, which stay lazy.
2.  Standardise it with
    [`el_decostand()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_decostand.html)
    (`total`, `hellinger`, `log`, `pa` or `rclr`), also lazy.
3.  Compute distances with
    [`el_vegdist()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_vegdist.html):
    `bray`, `jaccard`, `euclidean`, `manhattan`, `canberra`, `aitchison`
    or `robust.aitchison`. The result is a `dist`, and `vegan` takes
    over.
4.  Or compute the principal components with
    [`el_rda()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_rda.html),
    which returns the first `k` axes as an `rda` object.

``` r

# 1. ASVs present in at least 20% of the samples, and samples with at least 5,000 reads.
n_samples <- counts |> summarise(n = n_distinct(sample)) |> pull(n)

common <- counts |>
  summarise(prevalence = n_distinct(sample), .by = taxon) |>
  filter(prevalence >= 0.2 * n_samples) |>
  select(taxon)

deep <- counts |>
  summarise(depth = sum(quant, na.rm = TRUE), .by = sample) |>
  filter(depth >= 5000) |>
  select(sample)

filtered <- counts |>
  semi_join(common, by = "taxon") |>
  semi_join(deep, by = "sample")

# 2 and 3. Relative abundances and Bray-Curtis distances.
# Nothing is read into R until el_vegdist().
# `samples` sets the order of d to that of the metadata and includes samples without counts.
d <- filtered |>
  el_decostand("total") |>
  el_vegdist("bray", samples = metadata$sample)

# d is a dist object: use vegan on it.
ordination <- vegan::metaMDS(d)
vegan::adonis2(d ~ treatment, data = metadata)
vegan::envfit(ordination, metadata[, c("temperature", "salinity")], na.rm = TRUE)
```

`metaMDS()` on a `dist` has no species scores; giving it the matrix
would, but needs the dense table.

A PCA of the centred log-ratio (clr) of the counts needs a pseudocount,
because the logarithm of zero is not defined. `"inverse_depth"` uses one
over the sample’s total count, and a numeric vector named by sample sets
it freely. `"inverse_depth"` takes the totals of the table it is given,
which after the filters above are smaller than the libraries, by
different amounts in different samples. Pass the depths from before
filtering by name instead:

``` r

depth <- counts |>
  summarise(total = sum(quant, na.rm = TRUE), .by = sample) |>
  as.data.frame()
pseudocount <- purrr::set_names(1 / depth$total, depth$sample)

# 4. A PCA of the clr, first five axes.
pca <- filtered |>
  el_rda(k = 5, transform = "clr", pseudocount = pseudocount, samples = metadata$sample)

summary(pca)
plot(pca)

# Axis labels with the percentage of the total variance: divide by tot.chi, not by sum(eigenvals).
percent <- 100 * vegan::eigenvals(pca)[1:2] / pca$tot.chi
plot(
  vegan::scores(pca, display = "sites", choices = 1:2),
  xlab = sprintf("PC1 (%.1f%%)", percent[1]), ylab = sprintf("PC2 (%.1f%%)", percent[2])
)
vegan::envfit(pca, metadata[, c("temperature", "salinity")], choices = 1:2, na.rm = TRUE)
```

The prevalence filter matters for `clr`: on a 60,633-ASV time series
with 99% zeros, the seasonal cycle showed poorly in a `clr` PCA until
ASVs present in fewer than 20% of the samples were removed. `duckplyr`
translates few R functions, so `filter(.by = )` is not available (use
`summarise(.by = )` and a join, as above),
[`sum()`](https://rdrr.io/r/base/sum.html) needs `na.rm = TRUE`, and a
scalar with names, such as the result of
[`quantile()`](https://rdrr.io/r/stats/quantile.html), must go through
[`unname()`](https://rdrr.io/r/base/unname.html).

`clr` is not a method of
[`el_decostand()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_decostand.md),
because it makes every cell non-zero, which a table without zeros cannot
hold. It is computed from the raw counts inside the two functions that
use it: `el_vegdist(method = "aitchison")` is the Euclidean distance
between clr-transformed samples, which is the Aitchison distance, and
`el_rda(transform = "clr")` is a PCA of the clr, whose sample distances
with all axes are the same Aitchison distances.

A longer walk-through, from the published time series to seasonal
ordinations and alpha diversity, is in the vignette
[`seasonal-ordination`](https://biodiversitydata-se.github.io/ecolongr/articles/seasonal-ordination.html),
or
[`vignette("seasonal-ordination", package = "ecolongr")`](https://biodiversitydata-se.github.io/ecolongr/articles/seasonal-ordination.md).

## Where results differ from `vegan`

Distances, standardisations, diversity indices and expected richness
equal `vegan`’s to about 1e-10, and PCA scores to about 1e-8 (1e-6 with
the iterative solver on large tables). The tests compare every function
and option with `vegan` on `varespec`, `BCI` and a synthetic community.
The differences that matter:

- [`el_rrarefy()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_rrarefy.md)
  is random and does not give `vegan`’s subsample for the same seed; it
  draws from the same distribution.
- `el_decostand("rclr")` leaves zeros absent, which `vegan`’s completed
  `rclr` does not; `el_vegdist(method = "robust.aitchison")` and
  `el_rda(transform = "rclr")` run the same completion as `vegan`.
- `metaMDS()` on the `dist` has no species scores.
- [`el_rda()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_rda.md)
  returns the first `k` axes only, so
  [`predict()`](https://rdrr.io/r/stats/predict.html) and
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) give an
  error. Divide an eigenvalue by `pca$tot.chi`, never by the sum of the
  kept eigenvalues, to get the share of the whole variance.
- [`el_rarefy()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_drarefy.md)
  has no `se`, and a sample with no non-zero row is absent from the
  alpha diversity results, where `vegan` gives it a row of zeros.

All of them, with the reasons, are in the vignette
[`differences-from-vegan`](https://biodiversitydata-se.github.io/ecolongr/articles/differences-from-vegan.html).

## Memory and speed

[`el_vegdist()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_vegdist.md)
runs the lazy query in DuckDB, writes it to a temporary Parquet file
sorted by taxon, and reads it back in batches of whole taxa. A compiled
kernel adds up the pairs into one vector in the order of a `dist`
object. The memory has three parts:

- the result, one number per pair of samples, which grows with the
  square of the number of samples;
- one batch of rows;
- the memory DuckDB uses for the sort, which `memory_limit` caps (250 MB
  by default).

Measured on a 24-thread desktop under load, on synthetic tables where a
few taxa occur in nearly every sample:

| table | [`el_vegdist()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_vegdist.md) | [`el_rda()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_rda.md) |
|----|----|----|
| 5,000 samples, 11.1M rows | 7 s with 12 threads (30 s with one), 1.0 GB | 35 s, 1.4 GB |
| 20,000 samples, 19.4M rows | 25 s with 12 threads (148 s with one), 2.4 GB | 63 s, 2.2 GB |

Settings:

- `threads` (option `ecolongr.threads`, default half the cores) sets the
  threads of the distance kernel. The result is identical for any
  number.
- `memory_limit` (option `ecolongr.memory_limit`, default `"250MB"`)
  limits the memory DuckDB uses for sorting, and it spills to disk
  beyond that.
- `batch_rows` (default 1,000,000) is roughly how many table rows are in
  R at a time.

`ecolongr` does not change the settings of your `duckplyr` session. The
lazy steps you write before calling it, such as
[`el_decostand()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_decostand.md),
run there with DuckDB’s defaults (all cores, 80% of the memory). On a
small machine, limit them:

``` r

duckplyr::db_exec("SET memory_limit = '1GB'")
duckplyr::db_exec("SET threads = 4")
```

On the 5,000-sample table this took the peak of a `total` plus Euclidean
pipeline from 2.4 GB to 1.1 GB at a cost of 3 seconds.

## Development

- [`DESIGN.md`](https://biodiversitydata-se.github.io/ecolongr/DESIGN.md)
  has the scope, decisions, measurements and open questions.
- [`AGENTS.md`](https://biodiversitydata-se.github.io/ecolongr/AGENTS.md)
  has the working rules: tests first, `vegan` as the reference.
- `bench/` has scripts for the large synthetic tables. Benchmark with an
  installed package, not `devtools::load_all()`, which compiles without
  optimisation.
- `devtools::test()` runs the tests.

## License

MIT, see `LICENSE.md`.
