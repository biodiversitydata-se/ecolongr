# ecolongr (development version)

* Fixed `el_vegdist()` giving wrong distances, differing between runs, when `threads` was above 1 and the table had more than about 120,000 rows.
  DuckDB wrote the table sorted by taxon in parallel without keeping the order, so a taxon could be split in two and the pairs between its halves were lost.
  All methods were affected; `el_rda()` was not, because it does not depend on the order.
  The sorted table is now checked as it is read, so a wrong order is an error and not a wrong result.

# ecolongr 0.1.0

First version, for feedback.

* `el_vegdist()` computes `bray`, `jaccard`, `euclidean`, `manhattan`, `canberra`, `aitchison` and `robust.aitchison` dissimilarities from a long table without zeros, without building the taxa by samples matrix.
  The result is a `dist` and agrees with `vegan::vegdist()` to about 1e-10.
* `el_decostand()` makes `total`, `hellinger`, `log`, `pa` and `rclr` standardisations as a lazy long table.
  `clr` is computed from the counts inside `el_vegdist(method = "aitchison")` and `el_rda(transform = "clr")`, because it makes every cell non-zero.
* `el_rda()` returns an unconstrained PCA of the table, of its clr or of its completed robust clr, with the first `k` axes, as a `vegan` `rda` object.
  `scores()`, `eigenvals()`, `summary()`, `plot()`, `biplot()`, `envfit()` and `procrustes()` work on it, and `tot.chi` gives the variance the `k` axes leave out.
* `el_rda()` prints a message with the share of the total variance that the `k` axes hold and the share it leaves out.
  Silence it with `suppressMessages()` or `options(ecolongr.quiet = TRUE)`.
  The documentation explains why an eigenvalue is divided by `tot.chi` and not by the sum of the kept eigenvalues.
* The completion rank of the robust clr is `ropt`, 3 as in `vegan`, in `el_vegdist()` and `el_rda()`.
* `pseudocount` for `clr` is a number, `"inverse_depth"` (1 / the sample's total count) or a vector named by sample.
* A vignette follows the published LMO 16S time series from the Darwin Core archive to ordinations.
* Not implemented: constrained `rda` with explanatory variables, `predict()` and `fitted()` for `el_rda()` results, and species scores for `metaMDS()` on a `dist`.
* `R CMD check` runs in CI on Linux, macOS and Windows, with the tests and the vignette.
  The timings in the documentation were measured on Linux only.
