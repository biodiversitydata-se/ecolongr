# longvegan: design

Working name.
Status: draft for discussion, nothing implemented.

## Problem

`vegan` needs the full taxa × samples matrix in memory, dense.
For amplicon data with millions of ASVs and thousands of samples this is tens of gigabytes, although almost all cells are zero.
Selecting rows and columns first helps, but does not remove the problem.

Related: [asvoccur#1](https://github.com/biodiversitydata-se/asvoccur/issues/1) proposes keeping occurrence data as a long table (Parquet/DuckDB).
This package would compute the common `vegan` analyses directly from such a table.

## Goal

Re-implement a small subset of `vegan`, taking a long table (or a sparse matrix built from one) as input, with results that agree with `vegan` to numerical tolerance.
`vegan` is the reference implementation and the test oracle.

## Non-goals

- Replacing `vegan`.
  Anything that only needs a distance matrix is left to `vegan`.
- Constrained ordination, in the first versions.
- Matrix completion (`decostand(impute = TRUE)`).
- Plotting, including `ggplot2` interfaces.
  Result objects keep the `vegan` classes, so existing tools that plot `vegan` objects (`ggvegan`, `ggordiplots`, `microViz`) should work on them.
  Which of these are maintained and cover `rda` and `dist` has not been checked.

## Where the memory goes

The taxa axis is the large one.
After a sample × sample distance exists, that axis is gone.
So methods split by whether they need the taxa × samples matrix:

| vegan method | Needs the matrix? | Plan |
|---|---|---|
| `vegdist` | yes | implement |
| `decostand` | yes | implement the methods that need only non-zero cells |
| `rda`, unconstrained (PCA) | yes | implement with truncated SVD |
| `metaMDS` | no, takes a `dist` | use `vegan` on our `dist` |
| `adonis2` | no, takes a `dist` | use `vegan` on our `dist` |

`metaMDS` wraps `monoMDS`, a nonmetric optimiser with random starts and local minima, so a re-implementation would only agree approximately.
`adonis2` makes an n × n copy internally, which only matters at very large n.
Both can be revisited if profiling says otherwise.

## Input

The `duckplyr` direction is below; the exact input classes are open.
Candidates:

1. A long table with columns `sample`, `taxon`, `count`, as a data frame, an Arrow dataset or a DuckDB table.
2. A `dgCMatrix` (taxa × samples) built from (1) after selection.

A sparse matrix is the natural in-memory format for the linear algebra (`irlba`, `RSpectra`, C++ kernels).
The problem to avoid is the dense conversion, not the sparse matrix.
The intended flow is therefore: query the long table, select, build a sparse matrix only where an algorithm needs one.

## duckplyr interface

Users filter, join and reshape with `duckplyr` verbs, then pass the unevaluated table to `longvegan`.
Each function keeps working lazily for as long as the operation allows, and materializes only where it must.

### Contract

- **Input:** a lazy table (`duckplyr_df` or `duckdb_tibble`) with a sample, a taxon and a count column.
  The column names are arguments with defaults.
  A plain data frame is converted on entry.
- **`samples` (optional):** a lazy table or vector of sample IDs defining the full set of samples.
  A long table without zeros cannot show a sample whose counts are all zero, so such a sample is silently absent unless `samples` lists it.
  Without `samples`, the sample set is the samples that occur in the table, and the documentation says so.
  In `asvoccur` the `events` table is the natural source.
- **`decostand`:** lazy in, lazy out, as a long table.
- **`vegdist`:** the computation is lazy.
  The result is a `dist`, which is n × n and so materialized by definition.
- **`rda` (PCA):** pulls only the non-zero triplets into R to build the sparse matrix for `irlba`.
  The size is the number of non-zero cells, not taxa × samples.
- Each function documents where it materializes.

### Silent fallback

`duckplyr` falls back to `dplyr` without an error when it meets a verb it cannot translate, and the whole table is then materialized in R.
The only sign is the `duckplyr: materializing` message, shown with `options(duckdb.materialize_message = TRUE)`.
For a package whose purpose is saving memory this must not happen unnoticed, so:

- Every function needs a test that its pipeline stays lazy until the documented materialization point.
- A guard should turn an unexpected materialization into an error.
  How `duckplyr` allows this is not yet checked (see `?duckplyr::fallback_config`).

### Evidence so far

A scratch test on `vegan::varespec` as a long table:

- `semi_join`, `filter`, a many-to-many self-join on taxon, grouped `summarise`, further joins and `mutate` all stayed lazy until `as.data.frame()`.
- Bray–Curtis built from that pipeline matched `vegan::vegdist` to about 5e-16 over 66 sample pairs.
- `rclr` as `mutate(log(count) - mean(log(count)), .by = sample)` ran lazily.
  Its output has not yet been compared with `vegan`.
  The same formula in base R matches `decostand(..., "rclr", impute = FALSE)` exactly on the non-zero cells.
- `mutate(z = rev(count), .by = sample)` fell back to `dplyr` and materialized.

### Cost

The Bray–Curtis self-join produces, for each taxon, (number of samples containing it)² rows.
A taxon present in every sample costs n².
DuckDB aggregates before returning, so memory should be fine, but run time on large data is unmeasured.

## Function names

Recommendation: a prefix, for example `lv_vegdist`, `lv_decostand`, `lv_rda`.

- The functions are not drop-in replacements: the input is a lazy long table, only a subset of the options exists, and `rda` returns the first *k* axes.
- With the same names, `library(vegan); library(longvegan)` masks `vegdist`, `decostand` and `rda`.
  Old scripts then fail confusingly, or silently call the wrong function.
- `vegan::vegdist` is a plain function, not an S3 generic, so a method for our input class cannot coexist with it.
- Writing `longvegan::vegdist()` avoids the clash, but only if every example and user script does so.
- A prefix lets both packages be attached, and shows in each call which one is meant.
  It also groups the functions in autocompletion.
- The prefix should follow the package name, so it changes if the package is renamed.

Result objects keep the `vegan` classes (`dist`, and `rda` where feasible), so that `vegan::metaMDS`, `vegan::adonis2`, `scores()` and `plot()` work on them.

## Scope of the first version

### `vegdist`

- Indices: `bray`, `jaccard`, `euclidean`.
  Later: `manhattan`, `canberra`.
- Option: `binary`.
- Bray–Curtis is `1 - 2·Σ min(x_i, x_j) / (Σ x_i + Σ x_j)`.
  Only taxa present in both samples contribute to the sum of minima.
  Cost is the sum over taxa of (samples containing the taxon)².
- Returns a `dist` object, so the result can go straight into `vegan::metaMDS`, `vegan::adonis2`, `hclust` and so on.
- Samples with zero total counts are only seen if listed in `samples`.
  Listed ones are treated as all-zero rows, as in `vegan`:
  - Bray–Curtis and Jaccard are `NaN` for a zero–zero pair, and 1 for a zero sample against a non-empty one.
  - Euclidean between a zero sample and another is the norm of the other.
  - Omitting such a sample changes the result, not only the dimensions of the `dist`.

### `decostand`

Methods that act within a sample and need only the non-zero cells, so the result stays sparse:

- `total`, `hellinger`, `log`, `pa`.
- `rclr`: `log(x / g(x > 0))`, where `g` is the geometric mean of the positive values in the sample.
  `vegan` 2.7-6 returns `NA` for zeros with `impute = FALSE`, and the non-zero cells equal `log(x) - mean(log(x[x > 0]))` exactly.
  Called without `impute`, it fills the zeros by matrix completion, so the reference call must pass `impute = FALSE`.
  In a long table the zeros are simply absent rows.

Not planned: `clr` with a pseudocount, which fills the matrix.

### `rda`, unconstrained

- Truncated SVD of the centred (and optionally scaled) matrix, with the centring applied implicitly so the matrix is never densified (`irlba`, `RSpectra`).
- Total inertia is computed without densifying: `Σ x² − n · mean²` per column.
  This gives the proportion explained per axis.
- *n* is the number of samples in the full set.
  All-zero samples listed in `samples` count in *n*, the column means and the total inertia, as zero rows do in `vegan`.
- Only the first *k* axes are returned, where `vegan` returns all.
- To match: eigenvalues, site scores, species scores and `scaling` 1, 2 and 3.
- Axis signs are arbitrary.

## Testing against vegan

This is the main deliverable.

- Fixtures: `vegan::varespec`, `vegan::dune`, and a synthetic sparse community generated with a fixed seed.
- Compare with `testthat` at a tolerance of about 1e-10 for distances and transforms.
- Compare ordinations up to sign: eigenvalues, scores and proportion explained.
- Corner cases, each a separate test:
  - sample with all zeros, listed in `samples`: must match `vegan` on the matrix that includes the zero row
  - the same sample with `samples` omitted: must match `vegan` on the matrix without it
  - taxon present in a single sample
  - sample with a single taxon
  - samples with identical profiles
  - `binary = TRUE`
- Check that tests can fail: for each, a plausible wrong implementation (for example Bray–Curtis on relative abundances, or PCA without centring) must produce a failing test.
- Memory benchmark against `vegan` on a large synthetic dataset, reporting peak memory and run time.

## Existing packages to check first

Only the CRAN descriptions have been read so far, not the function lists.

- `proxyC`: proximity on large sparse matrices in C++.
  Not yet checked whether it has Bray–Curtis.
- `rbiom`: beta diversity, including UniFrac, on BIOM data.
- `irlba`, `RSpectra`: sparse truncated SVD.

If one of them already computes a measure correctly on sparse input, a thin wrapper plus the `vegan` comparison tests may be enough for that measure.

## Open questions

- Package name, and whether it belongs under the same organisation as `asvoccur`.
- Function names: prefix (recommended) or the `vegan` names.
  If a prefix, which one.
- Input contract: long table, sparse matrix, or both.
- Which lazy-table classes to accept: `duckplyr_df` and `duckdb_tibble` only, or also `tbl_lazy` from `dbplyr`, Arrow datasets and plain data frames.
- How to make an unexpected materialization an error, if `duckplyr` provides no switch for it.
- C++ (Rcpp) kernels, or DuckDB/`duckplyr` queries, for `vegdist`.
  A query is simpler to write, a kernel is probably faster.
- Whether to depend on `vegan`.
  It is needed for tests, and for returning compatible objects (`vegan::rda` class) if we choose to.
- Whether `rda` should return an object of class `rda`, so that `vegan`'s `scores`, `plot` and `summary` work.
- Whether omitting `samples` should only be documented, or also give a message that all-zero samples cannot be detected.
- Which indices beyond Bray–Curtis, Jaccard and Euclidean are worth the effort.
