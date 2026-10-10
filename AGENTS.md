# ecolongr

R package: community-ecology analyses (a subset of `vegan`) on long tables without zeros, computed lazily with `duckplyr`.
`DESIGN.md` holds the scope, the decisions and the open questions.
Read it first.

## Working rules

- Test first.
  Write the failing `testthat` test, show it, then implement.
- `vegan` is the reference: every result is compared with the `vegan` call it replaces, at tolerance 1e-10.
- Each test must be able to fail on a plausible wrong implementation.
- Functions use the prefix `el_` and take a table with columns `sample`, `taxon` and `quant`.

## duckplyr

- Convert input with `as_duckdb_tibble(x, prudence = "stingy")`, so an untranslatable verb is an error and not a silent copy into R.
  `as_duckdb_tibble()` on a `duckplyr_df` resets the prudence.
- Use `.by =`, not `group_by()`.
- Use `na.rm = TRUE` in `sum()` and `mean()` inside `mutate(.by = )`; without it `duckplyr` falls back.
- Collect only where `DESIGN.md` says a function materializes.
- Translatable: `log`, `log10`, `abs`, `round`, `as.integer`, arithmetic, comparisons, `if_else`, `min`, `max`, `sum`, `mean`.
  Not translatable: `sqrt`, `exp`, `floor`, `^`, `pmin`, `ifelse`, `cross_join`.
  Use `round(x) != x` for "not a whole number", a window `min()` in `mutate()` for a table-wide value, and `read_sql_duckdb()` after `compute_parquet()` where SQL is the only way.

## Commands

- `devtools::test()` runs the tests, `devtools::document()` regenerates `man/` and `NAMESPACE`, `devtools::check()` runs `R CMD check`.
- `duckplyr::fallback_sitrep()` shows the fallback logging settings.
  Never upload fallback reports.

## Performance and memory

- The kernel is `src/kernel.cpp` (Rcpp, OpenMP).
  Benchmark with an installed package: `devtools::load_all()` compiles with `-O0` and is about 10 times slower.
- `bench/gen.R` makes a skewed synthetic table and `bench/run.R` times `el_vegdist()` on it with peak memory.
- Run anything large on `sequoia` (24 cores, 125 GB, shared), not on the laptop (30 GB), and `nice` it.
  The repository is cloned at `~/dev/ecolongr` there.
- Peak memory is the reason for the design: the result (one number per pair of samples), one batch, and DuckDB's `memory_limit`.
  The query's first step runs in the `duckplyr` connection, whose settings belong to the user's session and must not be changed.
