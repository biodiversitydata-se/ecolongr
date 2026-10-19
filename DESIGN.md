# ecolongr: design

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

1. A long table with columns `sample`, `taxon`, `quant`, as a data frame, an Arrow dataset or a DuckDB table.
2. A `dgCMatrix` (taxa × samples) built from (1) after selection.

A sparse matrix is the natural in-memory format for the linear algebra (`irlba`, `RSpectra`, C++ kernels).
The problem to avoid is the dense conversion, not the sparse matrix.
The intended flow is therefore: query the long table, select, build a sparse matrix only where an algorithm needs one.

## Data sources

- [nf-core/ampliseq](https://nf-co.re/ampliseq) 3.0.0 and later write `summary_tables/ampliseq.counts.parquet`: a long table with columns `asv_id`, `sample` and `count`, zero rows dropped.
  It is DADA2 output before clustering and filtering, so it exists only with `--asv_calling dada2`.
- The same directory has one taxonomy table per classifier and database, keyed by `asv_id`, with flag columns such as `ampliseq_accept`.
  Join on `asv_id` to select ASVs, for example only those that passed all filters.
- `asvoccur` (see above) is the other source: public SBDI data, with `eventID` and `taxonID` in place of `sample` and `asv_id`.

Both give a long table with a sample, a taxon and a count, which is the input contract below.
A small subset of the ampliseq tables is a candidate test fixture beside the `vegan` data sets.

## duckplyr interface

Users filter, join and reshape with `duckplyr` verbs, then pass the unevaluated table to `ecolongr`.
Each function keeps working lazily for as long as the operation allows, and materializes only where it must.

### Contract

- **Input:** a lazy table (`duckplyr_df` or `duckdb_tibble`) with a sample, a taxon and a quantity column, by default `sample`, `taxon` and `quant`.
  The column names are arguments.
  The quantity is called `quant` and not `count`, because it may be a relative abundance or a transformed value.
  A plain data frame is converted on entry.
- **`samples` (optional):** a lazy table or vector of sample IDs defining the full set of samples.
  A long table without zeros cannot show a sample whose counts are all zero, so such a sample is silently absent unless `samples` lists it.
  Without `samples`, the sample set is the samples that occur in the table, and the documentation says so.
  In `asvoccur` the `events` table is the natural source.
- **`decostand`:** lazy in, lazy out, as a long table.
- **`vegdist`:** the lazy query is run by DuckDB and written to a temporary Parquet file.
  A compiled kernel reads it in batches of whole taxa and adds up the pairs.
  The result is a `dist`, which is n × n and so materialized by definition.
- **`rda` (PCA):** pulls only the non-zero triplets into R to build the sparse matrix for `irlba`.
  The size is the number of non-zero cells, not taxa × samples.
- Each function documents where it materializes.

### Example

```r
library(duckplyr)
library(ecolongr)

counts <- read_parquet_duckdb("summary_tables/ampliseq.counts.parquet")
taxonomy <- read_parquet_duckdb(
  "summary_tables/ampliseq.taxonomy.DADA2.consolidated_most-specific.parquet"
)

nmds <- counts |>
  semi_join(taxonomy, by = "asv_id") |>
  mutate(relab = count / sum(count, na.rm = TRUE), .by = sample) |>
  transmute(taxon = asv_id, sample, quant = relab) |>
  el_vegdist(method = "bray") |>
  vegan::metaMDS()
```

- The `semi_join()` keeps the ASVs that have a row in the taxonomy table.
  Taxonomy is assigned only after clustering and filtering, so this drops clustered-away and filtered ASVs from the counts.
  The file name is for `--consolidate_taxonomies most-specific`; with one database it is `ampliseq.taxonomy.DADA2.<database>.parquet`.
  A `filter()` on the taxonomy table, for example on `domain`, narrows the selection further before the join.
- Everything before `el_vegdist()` is lazy.
  `el_vegdist()` runs the query and returns a `dist`; `metaMDS()` is plain `vegan`.
- Use `.by =` for grouped work.
  `group_by()` is not translated by `duckplyr`: it falls back to `dplyr` and materializes.
- Window functions need `na.rm = TRUE`: `sum(count)` inside `mutate(.by = )` is not translated, `sum(count, na.rm = TRUE)` is.
  The same holds for `mean()`.
  The example above was checked on a stingy frame (see below) on 2026-10-10.
- The example renames to the default column names with `transmute()`.
- `ecolongr` does not attach `duckplyr`.
  Attaching `duckplyr` replaces `dplyr` verbs for data frames, which is a choice for the user's session, not for a package.

### Silent fallback

`duckplyr` falls back to `dplyr` without an error when it meets a verb it cannot translate, and the whole table is then materialized in R.
The only sign is the `duckplyr: materializing` message, shown with `options(duckdb.materialize_message = TRUE)`.
For a package whose purpose is saving memory this must not happen unnoticed, so:

- Every function needs a test that its pipeline stays lazy until the documented materialization point.
- `duckplyr` has no switch that turns a fallback into an error; `fallback_config()` only controls logging and uploading, by design.
  The guard is the `prudence` of the frame.
  A frame made with `as_duckdb_tibble(x, prudence = "stingy")` raises "This operation cannot be carried out by DuckDB" on any verb that would fall back, and `collect()` or `as.data.frame()` on a query DuckDB can run still work.
- Every function converts its input with `prudence = "stingy"`, so an untranslatable verb inside `ecolongr` is an error, not a silent copy.
- `as_duckdb_tibble()` on a frame that is already a `duckplyr_df` returns a lavish one, so the guard is applied by the function and does not depend on what the user passed.
- Users can do the same in their own pipelines by making the first table stingy.
  `ecolongr` cannot see a fallback that happened before the table reached it.

### Evidence so far

A scratch test on `vegan::varespec` as a long table:

- `semi_join`, `filter`, a many-to-many self-join on taxon, grouped `summarise`, further joins and `if_else` inside `summarise` run on a stingy frame.
- Bray–Curtis built from that pipeline matched `vegan::vegdist` to about 5e-16 over 66 sample pairs.
- `rclr` as `mutate(log(count) - mean(log(count), na.rm = TRUE), .by = sample)` runs on a stingy frame.
  Without `na.rm = TRUE` it falls back.
  An earlier version of this document said the version without `na.rm` ran lazily; that was judged from the result class and was wrong.
  Its output has not yet been compared with `vegan`.
  The same formula in base R matches `decostand(..., "rclr", impute = FALSE)` exactly on the non-zero cells.
- `mutate(relab = count / sum(count, na.rm = TRUE), .by = sample)` and `transmute()` run on a stingy frame; `group_by()` and `sum(count)` without `na.rm` do not.
- `mutate(z = rev(count), .by = sample)` fell back to `dplyr` and materialized.

### Cost

A many-to-many self-join on taxon produces, for each taxon, (number of samples containing it)² rows.
On a skewed synthetic table (a few taxa in nearly every sample) this did not work as a DuckDB query:

| samples | pair operations | DuckDB self-join | compiled kernel (first version) |
|---|---|---|---|
| 1,000 | 73M | 7.8 s, 2.9 GB | 0.7 s, 0.29 GB |
| 2,000 | 343M | killed at 14.6 GB after 30 min; 34 s and 4.3 GB with `memory_limit` | 3.9 s, 0.39 GB |

`el_vegdist()` therefore uses DuckDB only to select, sort and number the table, and a compiled kernel (`src/kernel.cpp`) to add up the pairs.
The kernel adds into one vector of n(n−1)/2 numbers, in the order of a `dist` object, so no pair table is built.
Samples are sorted within each taxon, blocks of sample ids are handled by one thread at a time, and the minima are computed in a vectorisable pass.
The result is bit-identical for any number of threads.

Measured on 5,000 samples, 11.1M table rows and 10.5 billion pair operations (a shared 24-core machine under load):

| version | time | peak memory |
|---|---|---|
| first kernel | 478 s | 1.8 GB |
| + sorting within taxa | 90 s | 1.6 GB |
| + branch-free minima, own DuckDB connection, integer ids (1 thread) | 30 s | 0.77 GB |
| same, 12 threads | 7.2 s | 1.0 GB |

The peak is the `dist` (100 MB here), one batch, DuckDB's `memory_limit` for the sort (default 250 MB, spills to disk), and the first step, which writes the lazy query to Parquet in the `duckplyr` connection (about 570 MB of buffers at 24 threads, which belongs to the user's session and is left alone).
Measured on 20,000 samples, 19.4M table rows and 60.5 billion pair operations: 25 s with 12 threads (148 s with one), peak 2.4 GB, of which the `dist` is 1.6 GB.
A dense taxa × samples matrix for the same table (500,000 taxa) would be about 80 GB.

Against `vegan::vegdist` on a smaller skewed table that does fit in memory (2,000 samples, 50,000 taxa, 928k rows, Bray–Curtis, from `bench/gen.R` with exponent 1), on a laptop (Intel Core i5-1245U, 12 threads, 30 GB):

| | time | peak memory |
|---|---|---|
| `el_vegdist()`, 4 threads | 1.0 s | 398 MB |
| `vegan::vegdist()` on the dense matrix (754 MB), 1 thread | not finished after about 14 minutes | 1.1 GB before it started |

`vegan::vegdist()` was stopped by a 15-minute limit on the whole script, so the real ratio is larger than 800.
Building the dense matrix from the long table is not counted in either time.
The two results were not compared on this table; equality with `vegan` is covered by the tests on smaller tables.
This is one skewed table and one metric, not a general speed-up.

Lazy steps before `el_vegdist()` run in the `duckplyr` connection of the user's session, with DuckDB's defaults (all cores, 80% of RAM).
On the 5,000-sample table, `el_decostand()` followed by `el_vegdist(method = "euclidean")` took 11–12 s and peaked at 2.4 GB for `total` and `rclr` and 3.9 GB for `log` (a table-wide window).
After `duckplyr::db_exec("SET memory_limit = '500MB'")` and `duckplyr::db_exec("SET threads = 4")` in the session, the same pipelines took 14–15 s and peaked at 1.1 GB.
The package does not change these settings; the documentation should tell users on small machines to set them.

Memory grows with the square of the number of samples: 1.6 GB for 20,000 samples and 10 GB for 50,000.

## Function names

Decision: the prefix `el_`, for example `el_vegdist`, `el_decostand`, `el_rda`.

- The functions are not drop-in replacements: the input is a lazy long table, only a subset of the options exists, and `rda` returns the first *k* axes.
- With the same names, `library(vegan); library(ecolongr)` masks `vegdist`, `decostand` and `rda`.
  Old scripts then fail confusingly, or silently call the wrong function.
- `vegan::vegdist` is a plain function, not an S3 generic, so a method for our input class cannot coexist with it.
- Writing `ecolongr::vegdist()` avoids the clash, but only if every example and user script does so.
- A prefix lets both packages be attached, and shows in each call which one is meant.
  It also groups the functions in autocompletion.

Result objects keep the `vegan` classes (`dist`, and `rda` where feasible), so that `vegan::metaMDS`, `vegan::adonis2`, `scores()` and `plot()` work on them.

## Scope of the first version

### `vegdist`

- Indices: `bray`, `jaccard`, `euclidean`, `manhattan`, `canberra` and `aitchison`.
  Manhattan is `Σx_i + Σx_j − 2 Σ min(x_i, x_j)` and uses the same kernel as Bray–Curtis.
  Canberra follows `vegan`: terms where both samples are zero are dropped and the sum is divided by the number of taxa in either sample, so it is `NA` between two all-zero samples.
  It needs the number of shared taxa as well, which the kernel counts in the same pass.
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

Methods that act within a sample and need only the non-zero cells, so the result stays a long table without zeros.
Implemented as `el_decostand()`; the result is a lazy `duckplyr` table and the tests compare every method with `vegan::decostand`.

- `total`, `log`, `pa`, `rclr`: lazy queries.
- `hellinger`: `duckplyr` cannot translate `sqrt()`, so the relative abundances are written to a temporary Parquet file and the square root is taken in SQL (`read_sql_duckdb()`).
  The file stays until the R session ends.
- `log` follows `vegan`: `log_b(x) + 1`, and if any value in the whole table is not a whole number, all values are first divided by the smallest positive value of the table (a table-wide window, not per sample).
  There is no warning.
- `rclr`: `log(x / g(x > 0))`, where `g` is the geometric mean of the positive values in the sample.
  `vegan` 2.7-6 returns `NA` for zeros with `impute = FALSE`, and the non-zero cells equal `log(x) - mean(log(x[x > 0]))` exactly.
  Called without `impute`, it fills the zeros by matrix completion, so the reference call must pass `impute = FALSE`.
  In a long table the zeros are simply absent rows.
- `el_decostand("rclr")` leaves zeros absent, so a later step that reads absent cells as 0 puts those taxa at the sample's geometric mean.
  That is not `vegan`'s robust Aitchison distance when zeros are present: `vegdist(method = "robust.aitchison")` is Euclidean distance on `decostand(x, "rclr")` with the default `impute = TRUE` (Martino et al. 2019).
  That step fills the missing cells by matrix completion of rank `ropt` (OptSpace, default rank 3, 5 iterations) and then centres rows and columns to mean 0, which also changes the non-zero cells.
  The two agree when the table has no zeros (tested).
  On the synthetic test table (80% zeros) they differ by up to 4.4 where the distances range up to 3.4.
- The completed version is implemented as `el_vegdist(method = "robust.aitchison")` and `el_rda(transform = "rclr")`, both from raw counts.
  OptSpace runs on the observed cells only (the algorithm, constants and line search follow `vegan::optspace()`, and the result agrees to 1e-14), with Rcpp helpers for the sums over cells.
  The completed matrix has rank `ropt`: after centring it is `Xc S Yc'`, so distances come from the `ropt`-column factors via a QR of `Yc`, and the PCA from an SVD of a `ropt` × `ropt` core.
  `el_rda` therefore returns at most `ropt` axes, `scale = TRUE` is not available, and samples with no non-zero cell are rejected (they have nothing to complete from).
- `ropt` is an argument of `el_vegdist()` (for `robust.aitchison`) and `el_rda()` (for `transform = "rclr"`), named and defaulted (3) as in `vegan`, where `decostand(x, "rclr", ropt = )` passes it to `optspace()`.
  It must be a whole number below the number of samples and of taxa, and giving it with any other method or transform is an error.
  Tested against `vegan` for ranks 2, 3 and 5.
  Speed and memory were measured only at rank 3: the cost of each OptSpace iteration grows with the rank, but how much on a large table is not measured.
  On the full LMO series (414 samples, 60,633 ASVs) it takes under two seconds, against 80 seconds for `vegan`'s dense completion, with the same seasonal R² (0.303).
  Measured on sequoia: 5,000 samples and 11.1M cells, `robust.aitchison` 38 s and 1.4 GB, `el_rda` 37 s and 1.7 GB; 20,000 samples and 19.4M cells, 87 s and 2.3 GB, and 72 s and 2.8 GB.
  The 20,000-sample run found an integer overflow in `p * n` (about 1e10), now fixed and covered by a test.
- `check = TRUE` (default) runs one query to stop on negative or missing quantities; `check = FALSE` is fully lazy.

`clr` with a pseudocount has no long-table form, because it makes every cell non-zero.
It is therefore not a method of `el_decostand()` (which stops and says so) but is computed from the raw counts inside `el_vegdist(method = "aitchison")` and `el_rda(transform = "clr")`.

### `clr`

- A clr row is `log(x + pc)` minus its mean over the *p* taxa.
  Write `b = log(1 + x / pc)`, which is 0 for a zero cell and sparse.
  The row mean of `log(x + pc)` is `log(pc)` plus the mean of `b`, so the clr row is `b` minus its row mean: a sparse matrix, centred per row.
- Aitchison distance (Euclidean on clr) between samples *i* and *j* is `sqrt(Σb_i² + Σb_j² − 2 Σ_shared b_i b_j − (Σb_i − Σb_j)² / p)`.
  The shared products are the same kernel as for Euclidean distance.
- The centred clr matrix is `B − u 1' − 1 c'` with `B` sparse, so `el_rda()` uses an implicit operator for `irlba`, and the exact solver for small tables.
- `pseudocount` is in the units of the counts: a positive number, `"inverse_depth"` for `1 / (the sample's total count)`, or a numeric vector named by sample.
  Give raw counts, not relative abundances.
  The user's practice is `count + 1/N_i` per sample, which compensates for the higher detection limit of shallow samples.
- The taxa of the table define *p*, as `vegan` after dropping all-zero columns.
  An all-zero sample listed in `samples` is the zero row of `vegan`: its clr is 0 and it sits at the origin.
- Measured on sequoia: 5,000 samples and 11.1M rows, `el_vegdist` 11 s and 1.0 GB, `el_rda` (5 axes) 35 s and 1.7 GB; 20,000 samples and 19.4M rows, 29 s and 2.8 GB, and 62 s and 3.1 GB.
- On the LMO time series (414 samples, 60,633 ASVs, 99% zeros), `el_rda` on the full table takes under a second and reproduces the dense `vegan` results exactly.

#### Choosing a transformation

On the LMO series, judged by how well PC1–PC2 show the seasonal cycle (see the memory notes of the analysis): Hellinger is good with or without filtering (seasonal R² 0.65), `clr` is poor on all 60,633 ASVs (0.28) but matches Hellinger once ASVs present in fewer than 20% of the samples are removed (0.64), and `rclr` with zeros as 0 is poor throughout (0.19–0.38).
The pseudocount mattered little, smaller ones doing slightly better.
`clr` should be used on a prevalence-filtered table.


`duckplyr` translates few functions.
Checked on 2026-10-10: `log`, `log10`, `abs`, `round`, `as.integer`, arithmetic and comparisons work; `sqrt`, `exp`, `floor`, `^`, `pmin`, `ifelse`, `cross_join` and `sum()` without `na.rm` inside a window do not.

### `rda`, unconstrained

Implemented as `el_rda()`.

- Truncated SVD of the centred (and optionally scaled) matrix, with the centring applied implicitly so the matrix is never densified.
  Small tables use a dense SVD (`solver = "exact"`), larger ones `irlba::irlba()` on a sparse matrix built from the non-zero rows.
- Total inertia is computed without densifying: `Σ x² − n · mean²` per column.
  This gives the proportion explained per axis.
- *n* is the number of samples in the full set.
  All-zero samples listed in `samples` count in *n*, the column means and the total inertia, as zero rows do in `vegan`.
- Only the first *k* axes are returned, where `vegan` returns all.
  The share of the total variance that they hold depends on *k*, so `el_rda()` reports it in a message after the fit, as `metaMDS()` reports its stress.
  `options(ecolongr.quiet = TRUE)` turns it off.
  The documentation says to divide an eigenvalue by `tot.chi`, not by the sum of the *k* eigenvalues, which sums to one and hides the rest.
- The result has class `ecolongr_rda`, `rda` and `cca`, so `scores()` (scalings 1, 2 and 3), `eigenvals()`, `summary()`, `plot()`, `biplot()`, `screeplot()`, `envfit()` and `procrustes()` work (tested against `vegan`).
  `predict()` and `fitted()` need the full matrix and give an error.
  `anova()` gives `vegan`'s table for an unconstrained ordination, which has no constrained component.
- The extra class exists for one `eigenvals()` method that sets the `sumev` attribute to the total inertia.
  Without it `summary()` takes the proportions of the *k* kept eigenvalues instead of the whole.
- Axis signs are arbitrary.
- Input is a transformed long table, composed with `el_decostand()`; `clr` with a pseudocount is not implemented.
  It fills the matrix, but only with one constant per sample, so it is a sparse matrix plus a rank-1 term and fits the implicit method.

Measured on sequoia, 5 axes, with 12 or fewer threads for the sort and `irlba` single-threaded:

| table | untransformed | `rclr` |
|---|---|---|
| 5,000 samples, 300,000 taxa, 11.1M rows | 35 s, 1.4 GB | 27 s, 3.0 GB |
| 20,000 samples, 499,000 taxa, 19.4M rows | 63 s, 2.2 GB | 54 s, 3.8 GB |

A dense matrix for the second table would be about 80 GB.
The higher peak with `rclr` is the lazy window step in the `duckplyr` connection (see above).

### `rda` with explanatory variables: postponed

`vegan::rda(Y ~ X)` takes a second, small matrix of explanatory variables, one row per sample, and splits the variation into a constrained and an unconstrained part.
It is not implemented, and the decision is to leave it for now.

What it would take:

- The maths fits the long-table approach.
  Regress the centred table on `X` through a QR of `X`: the fitted part has rank at most the number of variables, so its SVD is small, and the residual is the sparse table minus a low-rank part, which the implicit operator used for `clr` can handle with `irlba`.
- Partial RDA (`Condition()`) has the same structure.
- The cost is the permutation test.
  `vegan::anova()`, `permutest()`, `RsquareAdj()` and `ordiR2step()` read the dense fitted matrix stored in the object, which we do not have.
  We would write our own test: permuting the samples only needs the small matrix `Q'Y` recomputed, one sparse product per permutation, so it is feasible but is new code and a new object layout, not a flag on `el_rda()`.

What works now instead:

- `vegan::envfit()` fits variables to the axes of an unconstrained `el_rda()` result or of `metaMDS()` on an `el_vegdist()` distance, with permutation p-values (tested, including missing values with `na.rm = TRUE`).
  It is a post-hoc fit, not a constrained ordination.
- `vegan::dbrda()` and `capscale()` take the `dist` from `el_vegdist()`.
  The distance matrix grows with the square of the number of samples, which is fine up to a few thousand.
- `vegan::adonis2()` on a `dist` tests variables against the whole community.

Revisit it if constrained ordination of large tables, rather than of the distance matrix, turns out to be needed.

## Benchmark machine

Timings and peak memory in this document marked "measured on sequoia" come from one shared server.
They are indicative, not a specification: wall times depend on how busy the machine is, and peak memory does not.

| | |
|---|---|
| CPU | AMD Ryzen 9 3900X, 12 cores, 24 hardware threads, one socket |
| Memory | 125 GB RAM |
| Disk | NVMe SSD |
| System | Ubuntu 24.04, Linux 6.8 |
| Software | R 4.6.1, `duckplyr` 1.2.1, `duckdb` 1.5.4, `vegan` 2.7.5 |
| Load | Shared with other users' jobs; load average typically 16–24 on the 24 threads |

- Because the machine is busy, a run with `threads = 12` competes for cores and is slower than on an idle machine of the same size.
  Compare methods measured in the same session rather than across sessions.
- Peak memory is the resident set size (`VmHWM` in `/proc/self/status`) of the R process, which includes the DuckDB engine running inside it.
  Other users' jobs do not affect it.
- The package is installed with `R CMD INSTALL` and run with `Rscript`.
  `devtools::load_all()` compiles with `-O0` and is about ten times slower, so it is not used for timings.
- `bench/gen.R` writes a synthetic table as Parquet with a given number of samples, taxa and rows.
  Taxon frequencies follow a power law (Zipf), so a few taxa occur in nearly every sample.
  That is the expensive case for the pairwise kernel, whose work grows with the square of the number of samples per taxon.
  `bench/run.R` times `el_vegdist()` on such a file and reports the peak resident memory.
- The largest tables here (20,000 samples, 499,000 taxa) were run only on this server.
  A dense matrix of that table would take about 80 GB, which does not fit on a 30 GB laptop; the long-table methods needed 2–4 GB.

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

- Input contract: long table, sparse matrix, or both.
- Which lazy-table classes to accept: `duckplyr_df` and `duckdb_tibble` only, or also `tbl_lazy` from `dbplyr`, Arrow datasets and plain data frames.
- Whether to depend on `vegan`.
  It is needed for tests, and for returning compatible objects (`vegan::rda` class) if we choose to.
- Whether omitting `samples` should only be documented, or also give a message that all-zero samples cannot be detected.
- Whether to add `rda` with explanatory variables (see above).
- How the cost of the matrix completion grows with `ropt` above 3 on large tables, which has not been measured.
- Which indices beyond Bray–Curtis, Jaccard and Euclidean are worth the effort.
