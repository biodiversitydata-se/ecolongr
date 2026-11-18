# Alpha diversity of each sample

Counterparts of
[`vegan::diversity()`](https://vegandevs.github.io/vegan/reference/diversity.html)
and
[`vegan::specnumber()`](https://vegandevs.github.io/vegan/reference/diversity.html)
for a long table without zeros. The query runs in DuckDB and only the
result, one number per sample, is read into R; no taxa by samples matrix
is built. `el_diversity()` computes the sample totals, joins each total
back onto its rows and takes one grouped sum. `el_specnumber()` is one
grouped count.

## Usage

``` r
el_diversity(
  x,
  index = c("shannon", "simpson", "invsimpson"),
  base = exp(1),
  check = TRUE,
  sample = "sample",
  taxon = "taxon",
  quant = "quant"
)

el_specnumber(
  x,
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

- index:

  One of `"shannon"`, `"simpson"` and `"invsimpson"`.

- base:

  Base of the logarithm for `"shannon"`; ignored by the others.

- check:

  If `TRUE`, run one query to stop on negative, missing or non-finite
  (NaN, Inf) quantities. A sample total above about 1e308 is not
  checked. Use `FALSE` to skip the query.

- sample, taxon, quant:

  Names of the columns in `x`.

## Value

A numeric vector with one value per sample, named by sample and sorted
by name.

## Details

- `el_diversity()`: the Shannon index `-sum(p * log(p))` (divided by
  `log(base)`), the Simpson index `1 - sum(p * p)` or the inverse
  Simpson index `1 / sum(p * p)`, where `p` is the quantity divided by
  the sample total. As in `vegan`, quantities need not be whole numbers,
  and a common scale does not change the result, so proportions are
  accepted.

- `el_specnumber()`: the number of rows with a quantity above zero in
  each sample.

The table must have one row per sample and taxon. Rows are not merged,
so a sample and taxon that occur twice count as two taxa.

Rows with `quant == 0` are dropped. A sample with no non-zero row cannot
be told from a sample that is not in the table, so it is absent from the
result. `vegan` returns for an all-zero row Shannon 0, Simpson 1,
inverse Simpson `Inf` and richness 0; add these entries yourself if you
need them.

The names of the result are the samples as character, sorted as text:
factor samples are not in level order, numeric samples sort as
`"10" < "9"`, and an `NA` sample is named `NA` and comes last.

The sums are not compensated, so a sample with millions of taxa can
differ from `vegan` in the last digits (about 2e-11 relative with 1e7
equally common taxa), and which digits can vary between runs with
several threads.

Not implemented: `groups`, `MARGIN` and `equalize.groups` of
`diversity()`, and `renyi()`, `fisher.alpha()`, `estimateR()`,
`specpool()` and evenness.

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
richness <- el_specnumber(counts)
shannon <- el_diversity(counts, "shannon")
head(richness)
#> LMO_2011-2023_16S:P10353_1017 LMO_2011-2023_16S:P10353_1018 
#>                           406                           361 
#> LMO_2011-2023_16S:P10353_1019 LMO_2011-2023_16S:P10353_1020 
#>                           354                           433 
#> LMO_2011-2023_16S:P10353_1021 LMO_2011-2023_16S:P10353_1022 
#>                           313                           444 
head(shannon)
#> LMO_2011-2023_16S:P10353_1017 LMO_2011-2023_16S:P10353_1018 
#>                      4.161350                      3.411418 
#> LMO_2011-2023_16S:P10353_1019 LMO_2011-2023_16S:P10353_1020 
#>                      3.812296                      4.020182 
#> LMO_2011-2023_16S:P10353_1021 LMO_2011-2023_16S:P10353_1022 
#>                      3.428476                      4.273179 

# One table, with the metadata.
tibble(sample = names(shannon), richness = richness[names(shannon)], shannon = shannon) |>
  left_join(metadata, by = "sample") |>
  head()
#> # A tibble: 6 × 7
#>   sample                       richness shannon date  depth temperature salinity
#>   <chr>                           <int>   <dbl> <chr> <int>       <dbl>    <dbl>
#> 1 LMO_2011-2023_16S:P10353_10…      406    4.16 2017…     2        11.6      6.4
#> 2 LMO_2011-2023_16S:P10353_10…      361    3.41 2017…     2        13.6      6.6
#> 3 LMO_2011-2023_16S:P10353_10…      354    3.81 2017…     2        16.4      6.3
#> 4 LMO_2011-2023_16S:P10353_10…      433    4.02 2017…     2        14.6      6.2
#> 5 LMO_2011-2023_16S:P10353_10…      313    3.43 2017…     2        16.4      6  
#> 6 LMO_2011-2023_16S:P10353_10…      444    4.27 2017…     2        18.2      6.5

# Shannon index in bits.
el_diversity(counts, "shannon", base = 2) |> head()
#> LMO_2011-2023_16S:P10353_1017 LMO_2011-2023_16S:P10353_1018 
#>                      6.003559                      4.921636 
#> LMO_2011-2023_16S:P10353_1019 LMO_2011-2023_16S:P10353_1020 
#>                      5.499980                      5.799897 
#> LMO_2011-2023_16S:P10353_1021 LMO_2011-2023_16S:P10353_1022 
#>                      4.946245                      6.164894 
```
