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

## Commands

- `devtools::test()` runs the tests, `devtools::document()` regenerates `man/` and `NAMESPACE`, `devtools::check()` runs `R CMD check`.
- `duckplyr::fallback_sitrep()` shows the fallback logging settings.
  Never upload fallback reports.
