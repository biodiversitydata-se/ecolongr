utils::globalVariables("n_taxa")

#' Frequencies of frequencies
#'
#' For each sample, how many taxa were seen exactly once, twice, and so on.
#' This is the input of richness estimators that count the unseen taxa, such as
#' the `breakaway` package. It is one grouped count over the table, so it stays
#' lazy, and its size depends on the number of distinct counts, not the number
#' of taxa.
#'
#' @inheritParams el_rrarefy
#' @param check If `TRUE`, run queries that stop on negative, missing or
#'   non-whole counts. Use `FALSE` to stay entirely lazy.
#' @return A lazy `duckplyr` table with one row for each sample and count
#'   present in it: the sample column named as in `sample`, `count`, and
#'   `n_taxa`, the number of taxa with that count. A sample with no non-zero
#'   rows is absent. The row order is not defined.
#' @export
el_frequency_counts <- function(x, check = TRUE, sample = "sample", quant = "quant") {
  x <- duckplyr::as_duckdb_tibble(x, prudence = "stingy")
  x <- dplyr::select(x, sample = dplyr::all_of(sample), quant = dplyr::all_of(quant))
  if (check) check_counts(x, quant)
  x <- dplyr::filter(x, quant > 0)
  out <- dplyr::summarise(x, n_taxa = dplyr::n(), .by = c(sample, quant))
  out <- dplyr::select(out, sample, count = quant, n_taxa)
  dplyr::rename(out, dplyr::all_of(stats::setNames("sample", sample)))
}
