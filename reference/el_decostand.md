# Standardise a long table within samples

Lazy counterparts of the
[`vegan::decostand()`](https://vegandevs.github.io/vegan/reference/decostand.html)
methods that need only the non-zero cells. The result is a lazy
`duckplyr` table that can go straight into
[`el_vegdist()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_vegdist.md)
or further `duckplyr` verbs.

## Usage

``` r
el_decostand(
  x,
  method = c("total", "hellinger", "log", "pa", "rclr"),
  logbase = 2,
  check = TRUE,
  sample = "sample",
  taxon = "taxon",
  quant = "quant"
)
```

## Arguments

- x:

  A data frame or lazy `duckplyr` table with one row per non-zero sample
  and taxon. Rows with `quant == 0` are dropped.

- method:

  One of the methods above.

- logbase:

  Base of the logarithm for `method = "log"`.

- check:

  If `TRUE`, run one query to stop on negative, missing or non-finite
  (NaN, Inf) quantities. Use `FALSE` to stay entirely lazy.

- sample, taxon, quant:

  Names of the columns in `x`. The result has only these three columns,
  with the same names.

## Value

A lazy `duckplyr` table.

## Details

- `"total"`: quantities divided by the sample total.

- `"hellinger"`: square root of `"total"`.

- `"pa"`: 1 for every row.

- `"log"`: `log_b(x) + 1` as in `vegan`. If any value in the whole table
  is not a whole number, all values are first divided by the smallest
  positive value of the table, without a warning.

- `"rclr"`: `log(x)` minus the mean of `log(x)` over the sample's rows.
  This is `vegan`'s `"rclr"` with `impute = FALSE`; zeros stay absent
  where `vegan` returns `NA`. See the section on zeros below.

## Zeros and `rclr`

`vegan`'s default, `decostand(x, "rclr")`, fills the missing cells by
matrix completion and then centres rows and columns, which is a global
step that changes the non-zero cells too. `el_decostand()` does not:
zero cells stay absent, and a later step that reads absent cells as 0
(such as `el_vegdist(method = "euclidean")`) puts those taxa at the
sample's geometric mean. The result is therefore not `vegan`'s
`vegdist(method = "robust.aitchison")` when the table has zeros; the two
agree when it has none.

## Examples

``` r
library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union

# A long table of read counts for 16S ASVs in a time series from the Baltic Sea.
counts <- duckplyr::read_csv_duckdb(
  system.file("extdata", "lmo_free_living.tsv.gz", package = "ecolongr"),
  options = list(delim = "\t", header = TRUE), prudence = "stingy"
)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpqBfcfM/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
metadata <- read.delim(system.file("extdata", "lmo_free_living_samples.tsv", package = "ecolongr"))

# Samples with at least 5,000 reads.
deep <- counts |>
  summarise(reads = sum(quant, na.rm = TRUE), .by = sample) |>
  filter(reads >= 5000) |>
  select(sample)
counts <- semi_join(counts, deep, by = "sample")
metadata <- filter(metadata, sample %in% as.data.frame(deep)$sample)
# A lazy table: nothing is computed until a result is needed.
relative <- counts |> el_decostand("total")
relative |> head() |> as.data.frame()
#>                          sample                            taxon        quant
#> 1 LMO_2011-2023_16S:P10353_1017 00a01678e6e0adcd789980f8e123dfc5 3.427626e-05
#> 2 LMO_2011-2023_16S:P10353_1017 00db6fe233142d21d2f6c056132b95ed 9.793216e-06
#> 3 LMO_2011-2023_16S:P10353_1017 01c25b65bd637c1b1202bc70a73e6d8f 6.365591e-05
#> 4 LMO_2011-2023_16S:P10353_1017 01dff5d40f7a902fc54f6b7061676674 3.427626e-05
#> 5 LMO_2011-2023_16S:P10353_1017 02af743734992840c4f63515daf859ba 9.793216e-06
#> 6 LMO_2011-2023_16S:P10353_1017 02b019aa2218bf9b3d326950c6363127 1.542432e-03

# Standardise, then compute distances.
d <- counts |>
  el_decostand("hellinger") |>
  el_vegdist("euclidean", samples = metadata$sample)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpqBfcfM/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
attr(d, "Size")
#> [1] 335
```
