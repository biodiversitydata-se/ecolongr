# Where the results differ from vegan

Every function in `ecolongr` is compared with the `vegan` call it
replaces in the tests. Distances, standardisations, diversity indices
and expected richness equal `vegan`’s to about 1e-10, and PCA scores to
about 1e-8 (1e-6 with the iterative solver on large tables), apart from
the cases below. The tests use `varespec`, `BCI` and a synthetic
community, including all-zero samples and explicit zero rows.

## Alpha diversity

- **[`el_rrarefy()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_rrarefy.md)
  is random and does not give `vegan`’s subsample for the same seed.**
  [`vegan::rrarefy()`](https://vegandevs.github.io/vegan/reference/rarefy.html)
  draws the individuals one at a time.
  [`el_rrarefy()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_rrarefy.md)
  draws the count of each taxon from the hypergeometric distribution of
  what is left, which costs one draw per taxon. Both draw from the same
  distribution of counts, and the tests compare the mean and the
  variance of every cell over many replicates with `vegan`’s and with
  the exact values. The table is read into R, so this function needs a
  few times the memory of its input.
- **[`el_rarefy()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_drarefy.md)
  has no `se`, and takes one `depth`.** The standard error needs a sum
  over pairs of taxa. Call the function once for each depth.
- **A sample with no non-zero row is absent** from the result of
  [`el_specnumber()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_diversity.md),
  [`el_diversity()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_diversity.md)
  and
  [`el_rarefy()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_drarefy.md).
  `vegan` gives an all-zero row Shannon 0, Simpson 1, inverse Simpson
  `Inf` and richness 0.
- **`groups`, `MARGIN` and the other alpha-diversity functions**
  (`renyi()`, `fisher.alpha()`, `estimateR()`, `specpool()`) are not
  implemented.
- **The input has one row per sample and taxon.** A duplicated pair
  counts as two taxa, where `vegan` has already summed them in the
  matrix.

## Beta diversity

- **`el_decostand("rclr")` leaves zeros absent.** `vegan`’s
  `decostand(x, "rclr")` fills the zeros with values from a matrix
  completion, so every cell becomes non-zero. Like `clr`, that cannot be
  returned as a table without zeros. The completion is also one global
  fit over all observed cells, which a lazy table query cannot express.
  `el_decostand("rclr")` therefore does only the part that fits: the
  log-ratio of the non-zero cells. The completion also changes the
  non-zero cells; `el_decostand("rclr")` matches
  `decostand(x, "rclr", impute = FALSE)`. A later step that reads absent
  cells as 0 puts those taxa at the sample’s geometric mean, so
  Euclidean distance on it is not `vegan`’s `robust.aitchison` when the
  table has zeros. For that, use
  `el_vegdist(method = "robust.aitchison")` or
  `el_rda(transform = "rclr")` on the raw counts: they run the same
  matrix completion as `vegan`, of rank 3 unless you set `ropt`, on the
  observed cells only and agree with `vegan` to about 1e-6. They need
  every sample to have at least one non-zero cell.
- **`clr` needs a prevalence filter in practice.** On a 414-sample,
  60,633-ASV time series with 99% zeros, PC1–PC2 of a `clr` PCA showed
  the seasonal cycle poorly (seasonal R² 0.28), but as well as Hellinger
  once ASVs present in fewer than 20% of the samples were removed
  (0.64). Filter with `duckplyr` verbs before the call:
  [`semi_join()`](https://dplyr.tidyverse.org/reference/filter-joins.html)
  against the ASVs with
  `n_distinct(sample) >= 0.2 * n_distinct(all samples)`.
- **`metaMDS()` on the `dist` has no species scores.** Given a community
  matrix, `metaMDS()` also computes species scores as weighted averages
  and applies its own automatic transformations (`autotransform`). That
  needs the dense matrix, which is what this package avoids. On a `dist`
  you get the sample scores only: `ordination$species` is `NA` and
  `scores(ordination, display = "species")` is empty. Transform with
  [`el_decostand()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_decostand.md)
  before
  [`el_vegdist()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_vegdist.md)
  instead.
- **`log` follows `vegan` without the warning:** if any value is not a
  whole number, everything is first divided by the smallest positive
  value in the table.
- **`hellinger` writes a temporary Parquet file**, because `duckplyr`
  cannot compute a square root. The file is removed when the R session
  ends.
- **Bray–Curtis and Jaccard need non-negative quantities.** Euclidean
  distance and
  [`el_rda()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_rda.md)
  accept negative ones, for example after `rclr`.

## Ordination

- **[`el_rda()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_rda.md)
  returns the first `k` axes only.** `scores()`, `eigenvals()`,
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
  [`screeplot()`](https://rdrr.io/r/stats/screeplot.html) work.
  [`predict()`](https://rdrr.io/r/stats/predict.html) and
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) need the full
  matrix and give an error.
  [`anova()`](https://rdrr.io/r/stats/anova.html) gives `vegan`’s table
  for an unconstrained ordination, which has no constrained component.
  Axis signs are arbitrary, as in `vegan`.
- **`k` sets the number of axes (default 5).** A `k` above the rank,
  `min(number of samples − 1, number of taxa)`, is silently reduced to
  it, so a very large `k` returns all axes, as a full PCA. A `k` close
  to the rank makes
  [`el_rda()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_rda.md)
  use the exact solver, which builds the dense samples × taxa matrix.
  Keep `k` small on large tables, where the iterative solver is used.
  With `transform = "rclr"` there are at most `ropt` axes, the rank of
  the completion (default 3, as `vegan`’s `ropt`).
- **The variance the `k` axes leave out is known, and
  [`el_rda()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_rda.md)
  reports it in a message.** `pca$tot.chi` is the total variance of the
  whole table, whatever `k` is, so
  `1 - sum(vegan::eigenvals(pca)) / pca$tot.chi` is the unexplained
  share. Only that sum is known: the eigenvalues of the axes beyond `k`
  are not. Silence the message with
  [`suppressMessages()`](https://rdrr.io/r/base/message.html) or
  `options(ecolongr.quiet = TRUE)`.
- **Divide an eigenvalue by `tot.chi`, never by the sum of the
  eigenvalues.** `vegan::eigenvals(pca) / sum(vegan::eigenvals(pca))`
  sums to one over the `k` axes whatever `k` is, and so hides the
  variance that the other axes hold.
  `vegan::eigenvals(pca) / pca$tot.chi` and `summary(pca)` give the
  share of the whole.
- **Some plotting packages divide by the sum of the eigenvalues and so
  mislabel the axes.** `ggordiplots` (`ord_labels()`, used by its
  ordination plot functions) and `microViz` (`ord_plot()`) put a
  percentage in the axis label that is the eigenvalue divided by the sum
  of the kept eigenvalues, so the `k` labels sum to 100%. Labels from
  `eigenvals(pca) / pca$tot.chi`, as above, are right; replace the axis
  titles of such a plot with them, for example with `ggplot2::xlab()`.
  `ggvegan::autoplot()` does not print percentages and works on the
  result. For `transform = "rclr"` the share is 0, because the completed
  matrix has rank `ropt`, and the axes hold all of it.
