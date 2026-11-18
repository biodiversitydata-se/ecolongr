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
#' @param sample,quant Names of the columns in `x` that hold the sample and the count.
#' @param check If `TRUE`, run queries that stop on negative, missing or
#'   non-whole counts. Use `FALSE` to stay entirely lazy.
#' @return A lazy `duckplyr` table with one row for each sample and count
#'   present in it: the sample column named as in `sample`, `count`, and
#'   `n_taxa`, the number of taxa with that count. A sample with no non-zero
#'   rows is absent. The row order is not defined.
#' @example man-roxygen/lmo-counts.R
#' @examples
#' # How many ASVs each sample has seen once, twice and so on.
#' frequencies <- counts |> el_frequency_counts() |> as.data.frame()
#' head(frequencies)
#'
#' \dontrun{
#' # The input of richness estimators that count the unseen taxa, one sample at a time.
#' # `breakaway` is on GitHub (adw96/breakaway). It reads one table per sample, sorted by
#' # count, with the count in the first column and the number of taxa in the second.
#' per_sample <- frequencies |>
#'   arrange(sample, count) |>
#'   split(~ sample) |>
#'   purrr::map(\(d) select(d, count, n_taxa))
#'
#' estimates <- per_sample |>
#'   purrr::map(\(d) breakaway::breakaway(d))
#'
#' # The samples are independent, so the fits can run in parallel with purrr, which needs
#' # the mirai and carrier packages.
#' mirai::daemons(4)
#' estimates <- per_sample |>
#'   purrr::map(purrr::in_parallel(\(d) breakaway::breakaway(d)))
#' mirai::daemons(0)
#'
#' tibble(
#'   sample = names(estimates),
#'   estimate = purrr::map_dbl(estimates, \(e) e$estimate),
#'   se = purrr::map_dbl(estimates, \(e) e$error)
#' )
#' # breakaway's estimates depend strongly on `cutoff`, the largest count used in the fit.
#' # Tables without singletons, such as DADA2 output, need breakaway::breakaway_nof1().
#' }
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
