# Principal components of a long table

Unconstrained
[`vegan::rda()`](https://vegandevs.github.io/vegan/reference/cca.html)
on a table without zeros. The matrix is centred implicitly, so it is
never made dense, and only the first `k` axes are returned. Compose it
with
[`el_decostand()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_decostand.md)
for a transformation.

## Usage

``` r
el_rda(
  x,
  k = 5,
  scale = FALSE,
  samples = NULL,
  solver = c("auto", "exact", "irlba"),
  transform = c("none", "clr", "rclr"),
  pseudocount = 1,
  ropt = 3,
  sample = "sample",
  taxon = "taxon",
  quant = "quant"
)
```

## Arguments

- x:

  A data frame or lazy `duckplyr` table with one row per non-zero sample
  and taxon. Negative values are allowed.

- k:

  Number of axes to return, at most the rank, `min(n - 1, taxa)`. The
  axes beyond `k` are not computed, so the share of the total variance
  that the `k` axes hold depends on `k`. The function reports it in a
  message.

- scale:

  If `TRUE`, divide every taxon by its standard deviation, as in
  `vegan::rda(scale = TRUE)`.

- samples:

  Optional vector of all sample IDs. Samples in `samples` that are
  absent from `x` are all-zero rows and count in the means and the total
  inertia.

- solver:

  `"exact"` uses a dense SVD and is meant for small tables, `"irlba"` a
  sparse truncated SVD, and `"auto"` chooses by size.

- transform:

  `"none"` uses `x` as it is. `"clr"` takes raw counts and centred
  log-ratio transforms them, with every zero cell read as
  `log(pseudocount)`, without making the matrix dense. `"rclr"` takes
  raw counts, fills the zeros of the robust clr by matrix completion of
  rank `ropt` as in `vegan::decostand(x, "rclr")`, and returns at most
  `ropt` axes.

- pseudocount:

  For `transform = "clr"` only: a positive number, `"inverse_depth"` for
  `1 / (the sample's total count)`, or a positive numeric vector named
  by sample. The unit is that of the counts. `"inverse_depth"` uses the
  totals of the table passed in, so after filtering taxa it is not the
  depth of the whole library: compute the depths first and pass them as
  a named vector.

- ropt:

  For `transform = "rclr"` only: the rank of the matrix completion, 3 as
  in `vegan`. A whole number below the number of samples and of taxa.

- sample, taxon, quant:

  Names of the columns in `x`.

## Value

An object of class `rda` and `cca` with the first `k` axes. Use
[`vegan::scores()`](https://vegandevs.github.io/vegan/reference/scores.html)
and
[`vegan::eigenvals()`](https://vegandevs.github.io/vegan/reference/eigenvals.html)
on it. Functions that need the full fitted matrix, such as
[`predict()`](https://rdrr.io/r/stats/predict.html) and
[`anova()`](https://rdrr.io/r/stats/anova.html), do not work.

`tot.chi` is the total variance of the whole table, whatever `k` is. The
share of the variance that an axis holds is its eigenvalue divided by
`tot.chi`, which [`summary()`](https://rdrr.io/r/base/summary.html)
does. Dividing by the sum of the `k` kept eigenvalues is wrong: those
shares always sum to one and hide the variance that the other axes hold.

A message gives the share held by the `k` axes and the share left out.
[`suppressMessages()`](https://rdrr.io/r/base/message.html) or
`options(ecolongr.quiet = TRUE)` silences it.

## Examples

``` r
data(varespec, package = "vegan")
x <- data.frame(
  sample = rep(rownames(varespec), ncol(varespec)),
  taxon = rep(colnames(varespec), each = nrow(varespec)),
  quant = unlist(varespec, use.names = FALSE)
)
x <- x[x$quant > 0, ]
pca <- el_rda(x, k = 3)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/Rtmp5qBYmR/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> `el_rda()`: the first 3 axes hold 86.5% of the total variance and 13.5% is left out. Use `summary()` for the share of each axis.

# Right: the share of the whole variance. The three do not sum to one.
eig <- vegan::eigenvals(pca)
c(eig / pca$tot.chi)
#>        PC1        PC2        PC3 
#> 0.53842399 0.25432128 0.07243987 
summary(pca)$cont$importance["Proportion Explained", ]
#>        PC1        PC2        PC3 
#> 0.53842399 0.25432128 0.07243987 
1 - sum(eig) / pca$tot.chi
#> [1] 0.1348149

# Wrong: the three always sum to one, whatever the axes leave out.
c(eig / sum(eig))
#>        PC1        PC2        PC3 
#> 0.62232229 0.29395012 0.08372759 
```
