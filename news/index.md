# Changelog

## ecolongr (development version)

- Every function has examples in its help page, run on the table bundled
  in the package. The README is organised by application area (alpha
  diversity, beta diversity and ordination), the comparison with `mia`
  is removed, and the list of differences from `vegan` is a new
  vignette. There is a package website, built by `pkgdown` from the
  `main` branch.

- New
  [`el_frequency_counts()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_frequency_counts.md):
  for each sample, how many taxa were seen exactly once, twice and so
  on. It is the input of richness estimators that count the unseen taxa,
  such as the `breakaway` package, and stays lazy.

- New
  [`el_diversity()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_diversity.md)
  and
  [`el_specnumber()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_diversity.md),
  the counterparts of
  [`vegan::diversity()`](https://vegandevs.github.io/vegan/reference/diversity.html)
  (Shannon, Simpson and inverse Simpson, with `base`) and
  [`vegan::specnumber()`](https://vegandevs.github.io/vegan/reference/diversity.html):
  one value per sample, as a named vector. They are lazy, and only the
  per-sample result is read into R. On a table of 10.9M rows, 2,000
  samples and 30,000 taxa,
  [`el_diversity()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_diversity.md)
  added 37 MB of peak memory with 4 threads and 217 MB with 12, in
  0.2-0.3 s. A sample with no non-zero row is absent from the result,
  where `vegan` gives an all-zero row Shannon 0, Simpson 1, inverse
  Simpson `Inf` and richness 0. `groups`, `MARGIN`, `equalize.groups`
  and the other alpha-diversity functions are not implemented.

- The `check` of
  [`el_decostand()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_decostand.md),
  [`el_drarefy()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_drarefy.md),
  [`el_rarefy()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_drarefy.md),
  [`el_diversity()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_diversity.md)
  and
  [`el_specnumber()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_diversity.md)
  now also stops on `NaN` and `Inf` quantities, which it let through
  before.

- New
  [`el_drarefy()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_drarefy.md)
  and
  [`el_rarefy()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_drarefy.md),
  the counterparts of
  [`vegan::drarefy()`](https://vegandevs.github.io/vegan/reference/rarefy.html)
  and
  [`vegan::rarefy()`](https://vegandevs.github.io/vegan/reference/rarefy.html):
  the expected presence of each taxon, and the expected richness of each
  sample, in a subsample of a common depth. They equal `vegan`’s results
  and need no random numbers.
  [`el_rarefy()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_drarefy.md)
  has no `se` and takes one `depth`.

- New
  [`el_rrarefy()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_rrarefy.md),
  the counterpart of
  [`vegan::rrarefy()`](https://vegandevs.github.io/vegan/reference/rarefy.html):
  a random subsample of each sample to a common depth. It draws from the
  same distribution as `vegan` but not the same random numbers, so a
  seed does not give `vegan`’s result; the tests compare means and
  variances over many replicates with `vegan` and with the exact values.
  The table is read into R.

## ecolongr 0.1.1

- Fixed
  [`el_vegdist()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_vegdist.md)
  giving wrong distances, differing between runs, when `threads` was
  above 1 and the table had more than about 120,000 rows. DuckDB wrote
  the table sorted by taxon in parallel without keeping the order, so a
  taxon could be split in two and the pairs between its halves were
  lost. All methods were affected;
  [`el_rda()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_rda.md)
  was not, because it does not depend on the order. The sorted table is
  now checked as it is read, so a wrong order is an error and not a
  wrong result.
- The README compares the time and memory of
  [`el_vegdist()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_vegdist.md)
  with `vegan` and `mia`, and `DESIGN.md` is brought up to date.

## ecolongr 0.1.0

First version, for feedback.

- [`el_vegdist()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_vegdist.md)
  computes `bray`, `jaccard`, `euclidean`, `manhattan`, `canberra`,
  `aitchison` and `robust.aitchison` dissimilarities from a long table
  without zeros, without building the taxa by samples matrix. The result
  is a `dist` and agrees with
  [`vegan::vegdist()`](https://vegandevs.github.io/vegan/reference/vegdist.html)
  to about 1e-10.
- [`el_decostand()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_decostand.md)
  makes `total`, `hellinger`, `log`, `pa` and `rclr` standardisations as
  a lazy long table. `clr` is computed from the counts inside
  `el_vegdist(method = "aitchison")` and `el_rda(transform = "clr")`,
  because it makes every cell non-zero.
- [`el_rda()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_rda.md)
  returns an unconstrained PCA of the table, of its clr or of its
  completed robust clr, with the first `k` axes, as a `vegan` `rda`
  object. `scores()`, `eigenvals()`,
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html),
  [`biplot()`](https://rdrr.io/r/stats/biplot.html), `envfit()` and
  `procrustes()` work on it, and `tot.chi` gives the variance the `k`
  axes leave out.
- [`el_rda()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_rda.md)
  prints a message with the share of the total variance that the `k`
  axes hold and the share it leaves out. Silence it with
  [`suppressMessages()`](https://rdrr.io/r/base/message.html) or
  `options(ecolongr.quiet = TRUE)`. The documentation explains why an
  eigenvalue is divided by `tot.chi` and not by the sum of the kept
  eigenvalues.
- The completion rank of the robust clr is `ropt`, 3 as in `vegan`, in
  [`el_vegdist()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_vegdist.md)
  and
  [`el_rda()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_rda.md).
- `pseudocount` for `clr` is a number, `"inverse_depth"` (1 / the
  sample’s total count) or a vector named by sample.
- A vignette follows the published LMO 16S time series from the Darwin
  Core archive to ordinations.
- Not implemented: constrained `rda` with explanatory variables,
  [`predict()`](https://rdrr.io/r/stats/predict.html) and
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) for
  [`el_rda()`](https://biodiversitydata-se.github.io/ecolongr/reference/el_rda.md)
  results, and species scores for `metaMDS()` on a `dist`.
- `R CMD check` runs in CI on Linux, macOS and Windows, with the tests
  and the vignette. The timings in the documentation were measured on
  Linux only.
