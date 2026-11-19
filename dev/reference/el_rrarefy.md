# Rarefy samples to a common depth

Counterpart of
[`vegan::rrarefy()`](https://vegandevs.github.io/vegan/reference/rarefy.html):
a random subsample, without replacement, of `depth` individuals from
each sample. With no `depth`, the smallest sample total is used, so no
sample is dropped.

## Usage

``` r
el_rrarefy(
  x,
  depth = NULL,
  sample = "sample",
  taxon = "taxon",
  quant = "quant"
)
```

## Arguments

- x:

  A data frame or lazy `duckplyr` table with one row per non-zero sample
  and taxon, holding whole counts. Rows with `quant == 0` are dropped.

- depth:

  A single finite whole number, at least 1. The default, `NULL`, is the
  smallest total among the samples with a non-zero count.

- sample, taxon, quant:

  Names of the columns in `x`. The result has only these three columns,
  with the same names.

## Value

A `duckplyr` table holding the result in memory, with the samples in the
order they first appear in `x` and without the taxa that were not drawn.

## Details

A sample with a total below `depth` cannot be rarefied and is dropped
from the result, with a warning that names it; it is an error if no
sample is left. `vegan` returns such a sample unchanged, which leaves
samples of different depths in the result. A sample whose total equals
`depth` is kept as it is.

The result is random, so call
[`set.seed()`](https://rdrr.io/r/base/Random.html) first for a
reproducible one. It is not the same random sample as `vegan`'s for the
same seed: `vegan` draws the individuals one at a time, and
`el_rrarefy()` draws the count of each taxon from the hypergeometric
distribution of what is left. Both draw from the same distribution of
counts (the multivariate hypergeometric).

The table is read into R and copied several times, so this is not a lazy
step: allow a few times the memory of the input. Counts or totals of
about 2^31 and above make each draw slow.

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
set.seed(1)
subsample <- el_rrarefy(counts, depth = 5000)

# Every sample now has exactly 5,000 reads, but different samples keep different ASVs.
subsample |> summarise(reads = sum(quant), .by = sample) |> head() |> as.data.frame()
#>                          sample reads
#> 1 LMO_2011-2023_16S:P10353_1023  5000
#> 2 LMO_2011-2023_16S:P10353_1026  5000
#> 3  LMO_2011-2023_16S:P1202_1051  5000
#> 4 LMO_2011-2023_16S:P13006_1026  5000
#> 5 LMO_2011-2023_16S:P13006_1029  5000
#> 6 LMO_2011-2023_16S:P13006_1031  5000
el_specnumber(subsample) |> head()
#> LMO_2011-2023_16S:P10353_1017 LMO_2011-2023_16S:P10353_1018 
#>                           251                           236 
#> LMO_2011-2023_16S:P10353_1019 LMO_2011-2023_16S:P10353_1020 
#>                           222                           301 
#> LMO_2011-2023_16S:P10353_1021 LMO_2011-2023_16S:P10353_1022 
#>                           221                           303 
```
