# Dissimilarity from a long table

Computes the distance from a table without zeros. The query is run by
DuckDB and written to a temporary Parquet file sorted by taxon. A
compiled kernel then reads it in batches and adds up the shared
quantities for every pair of samples, so the memory used is the result
(one number per pair of samples) plus one batch.

## Usage

``` r
el_vegdist(
  x,
  method = c("bray", "jaccard", "euclidean", "aitchison", "manhattan", "canberra",
    "robust.aitchison"),
  binary = FALSE,
  pseudocount = 1,
  ropt = 3,
  samples = NULL,
  batch_rows = 1e+06,
  threads = NULL,
  memory_limit = NULL,
  sample = "sample",
  taxon = "taxon",
  quant = "quant"
)
```

## Arguments

- x:

  A data frame or lazy `duckplyr` table with one row per non-zero sample
  and taxon. Rows with `quant == 0` are ignored, but a sample that has
  only such rows is kept as an all-zero sample. Negative values are
  allowed only for `method = "euclidean"` without `binary`, for example
  after `el_decostand(method = "rclr")`.

- method:

  One of `"bray"`, `"jaccard"`, `"euclidean"`, `"aitchison"`,
  `"manhattan"` and `"canberra"`, as in
  [`vegan::vegdist()`](https://vegandevs.github.io/vegan/reference/vegdist.html).
  `"aitchison"` is the Euclidean distance between clr-transformed
  samples, computed from the counts. `"canberra"` is `NA`, as in
  `vegan`, between two all-zero samples. `"robust.aitchison"` is the
  Euclidean distance on the robust clr with the zeros filled by matrix
  completion (rank `ropt`, as in `vegan`). It needs counts, and no
  all-zero samples.

- binary:

  If `TRUE`, use presence and absence instead of quantities.

- pseudocount:

  For `method = "aitchison"` only: added to the counts before the clr. A
  positive number, `"inverse_depth"` for
  `1 / (the sample's total count)`, or a positive numeric vector named
  by sample. The unit is that of the counts, so give raw counts, not
  relative abundances. `"inverse_depth"` uses the totals of the table
  passed in, so after filtering taxa it is not the depth of the whole
  library: compute the depths first and pass them as a named vector.

- ropt:

  For `method = "robust.aitchison"` only: the rank of the matrix
  completion, 3 as in `vegan`. A whole number below the number of
  samples and of taxa.

- samples:

  Optional vector of all sample IDs, which also sets the order of the
  result. Samples in `samples` that are absent from `x` are treated as
  all zero. Without it, the samples that occur in `x` are used, in
  alphabetical order.

- batch_rows:

  Approximate number of table rows held in R at a time. A batch always
  holds whole taxa.

- threads:

  Number of threads for the pair kernel. Defaults to the option
  `ecolongr.threads`, or half of the cores. The result does not depend
  on it.

- memory_limit:

  Memory DuckDB may use for sorting the table, such as `"250MB"`; it
  spills to disk beyond that. Defaults to the option
  `ecolongr.memory_limit`, or `"250MB"`. This is separate from the
  memory used by R for the result and the current batch.

- sample, taxon, quant:

  Names of the columns in `x`.

## Value

A `dist` object.

## Details

Euclidean distance on `el_decostand(method = "rclr")` output reads
absent cells as 0. That is not `vegan`'s `"robust.aitchison"` when the
table has zeros, because `vegan` first fills them by matrix completion;
see
[`el_decostand()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_decostand.md).

## Examples

``` r
library(dplyr)

# A long table of read counts for 16S ASVs in a time series from the Baltic Sea.
counts <- duckplyr::read_csv_duckdb(
  system.file("extdata", "lmo_free_living.tsv.gz", package = "ecolongr"),
  options = list(delim = "\t", header = TRUE), prudence = "stingy"
)
metadata <- read.delim(system.file("extdata", "lmo_free_living_samples.tsv", package = "ecolongr"))

# Samples with at least 5,000 reads.
deep <- counts |>
  summarise(reads = sum(quant, na.rm = TRUE), .by = sample) |>
  filter(reads >= 5000) |>
  select(sample)
counts <- semi_join(counts, deep, by = "sample")
metadata <- filter(metadata, sample %in% as.data.frame(deep)$sample)
d <- counts |>
  el_decostand("total") |>
  el_vegdist("bray", samples = metadata$sample)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpXIBetp/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.

# A dist object, so base R and vegan work on it.
hc <- hclust(d, "average")
head(cutree(hc, k = 4))
#>  LMO_2011-2023_16S:P2256_1001 LMO_2011-2023_16S:P17302_1014 
#>                             1                             1 
#>  LMO_2011-2023_16S:P2256_1002 LMO_2011-2023_16S:P17302_1015 
#>                             1                             1 
#> LMO_2011-2023_16S:P17302_1016  LMO_2011-2023_16S:P2256_1003 
#>                             1                             1 
head(cmdscale(d, k = 2))
#>                                      [,1]        [,2]
#> LMO_2011-2023_16S:P2256_1001  -0.36373842  0.09539471
#> LMO_2011-2023_16S:P17302_1014 -0.25660013  0.09419716
#> LMO_2011-2023_16S:P2256_1002  -0.38406749  0.01317695
#> LMO_2011-2023_16S:P17302_1015  0.04824488 -0.23368001
#> LMO_2011-2023_16S:P17302_1016 -0.24735047  0.05039300
#> LMO_2011-2023_16S:P2256_1003  -0.38311115 -0.06718847

# Aitchison distance from the raw counts, with a pseudocount of 1 / library size.
d_aitchison <- counts |>
  el_vegdist("aitchison", pseudocount = "inverse_depth", samples = metadata$sample)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpXIBetp/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
round(cor(d, d_aitchison), 2)
#> [1] 0.35
```
