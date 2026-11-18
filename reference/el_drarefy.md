# Expected presence and richness at a common depth

Counterparts of
[`vegan::drarefy()`](https://vegandevs.github.io/vegan/reference/rarefy.html)
and
[`vegan::rarefy()`](https://vegandevs.github.io/vegan/reference/rarefy.html):
what a random subsample of `depth` individuals from each sample is
expected to contain. There is no randomness, and the results equal
`vegan`'s.

## Usage

``` r
el_drarefy(
  x,
  depth,
  check = TRUE,
  sample = "sample",
  taxon = "taxon",
  quant = "quant"
)

el_rarefy(
  x,
  depth,
  check = TRUE,
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

  A single finite whole number, at least 1.

- check:

  If `TRUE`, run queries that stop on negative, missing, non-finite or
  non-whole counts and warn about samples below `depth`. Use `FALSE` to
  stay entirely lazy.

- sample, taxon, quant:

  Names of the columns in `x`. The result has only these three columns,
  with the same names.

## Value

`el_drarefy()`: a lazy `duckplyr` table with the columns `sample`,
`taxon` and `quant`, named as in the arguments. `el_rarefy()`: a numeric
vector with one expected richness per sample, named by sample.

## Details

- `el_drarefy()`: for each non-zero cell, the probability that the taxon
  is in the subsample, `1 - choose(J - x, depth) / choose(J, depth)`
  with `J` the sample total. Absent taxa have probability 0 and stay
  absent, so the result is a lazy long table without zeros.

- `el_rarefy()`: the sum of these probabilities over a sample, which is
  its expected richness at `depth`. Unlike `vegan`, it has no `se`
  option and takes a single `depth`.

A sample with a total below `depth` is treated as if `depth` were its
total, so every taxon in it has probability 1, with a warning as in
`vegan`.

The probability depends only on the sample total and the count, so it is
computed once for each distinct pair of the two, in R with
[`lchoose()`](https://rdrr.io/r/base/Special.html) as `vegan` does. Only
these pairs, not the table, are read into R.

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
# Expected richness of each sample in a subsample of 5,000 reads.
expected <- el_rarefy(counts, depth = 5000)
observed <- el_specnumber(counts)
head(cbind(observed, expected = expected[names(observed)]))
#>                               observed expected
#> LMO_2011-2023_16S:P10353_1017      406 255.0769
#> LMO_2011-2023_16S:P10353_1018      361 235.8478
#> LMO_2011-2023_16S:P10353_1019      354 234.3010
#> LMO_2011-2023_16S:P10353_1020      433 304.5102
#> LMO_2011-2023_16S:P10353_1021      313 226.4471
#> LMO_2011-2023_16S:P10353_1022      444 296.1962

# The probability that each ASV is in such a subsample, as a lazy table.
counts |>
  el_drarefy(depth = 5000) |>
  head() |>
  as.data.frame()
#>                          sample                            taxon     quant
#> 1 LMO_2011-2023_16S:P10353_1017 272f79649c2702bec9ff272552ffffcf 0.7104818
#> 2 LMO_2011-2023_16S:P10353_1017 96acf2340c6d1dfb87344ddd9bcc9318 0.8403265
#> 3 LMO_2011-2023_16S:P10353_1017 f5245453d01e9ff275d2d3e0b2f3f270 0.9295593
#> 4 LMO_2011-2023_16S:P10353_1018 f3931878c88c9d308370e5e20aa824b3 0.3142853
#> 5 LMO_2011-2023_16S:P10353_1018 2a1ad2454c36352c90d72b9f35dc37d0 0.4929517
#> 6 LMO_2011-2023_16S:P10353_1018 c1363f1fd0e8a3fee5eecfe66c0ae7dc 0.9701010
```
