utils::globalVariables(c("total", "k", "p"))

#' Expected presence and richness at a common depth
#'
#' Counterparts of [vegan::drarefy()] and [vegan::rarefy()]: what a random
#' subsample of `depth` individuals from each sample is expected to contain.
#' There is no randomness, and the results equal `vegan`'s.
#'
#' * `el_drarefy()`: for each non-zero cell, the probability that the taxon is
#'   in the subsample, `1 - choose(J - x, depth) / choose(J, depth)` with `J`
#'   the sample total. Absent taxa have probability 0 and stay absent, so the
#'   result is a lazy long table without zeros.
#' * `el_rarefy()`: the sum of these probabilities over a sample, which is its
#'   expected richness at `depth`. Unlike `vegan`, it has no `se` option and
#'   takes a single `depth`.
#'
#' A sample with a total below `depth` is treated as if `depth` were its total,
#' so every taxon in it has probability 1, with a warning as in `vegan`.
#'
#' The probability depends only on the sample total and the count, so it is
#' computed once for each distinct pair of the two, in R with [lchoose()] as
#' `vegan` does. Only these pairs, not the table, are read into R.
#'
#' @inheritParams el_rrarefy
#' @param check If `TRUE`, run queries that stop on negative, missing, non-finite
#'   or non-whole counts and warn about samples below `depth`. Use `FALSE` to stay
#'   entirely lazy.
#' @return `el_drarefy()`: a lazy `duckplyr` table with the columns `sample`,
#'   `taxon` and `quant`, named as in the arguments. `el_rarefy()`: a numeric
#'   vector with one expected richness per sample, named by sample.
#' @export
el_drarefy <- function(x, depth, check = TRUE, sample = "sample", taxon = "taxon", quant = "quant") {
  x <- rarefy_input(x, depth, check, sample, taxon, quant)
  pairs <- as.data.frame(dplyr::distinct(x, total, quant))
  pairs$p <- presence_probability(pairs$total, pairs$quant, depth)
  pairs <- duckplyr::as_duckdb_tibble(pairs, prudence = "stingy")
  out <- dplyr::inner_join(x, pairs, by = c("total", "quant"))
  out <- dplyr::select(out, sample, taxon, quant = p)
  dplyr::rename(out, dplyr::all_of(stats::setNames(c("sample", "taxon", "quant"), c(sample, taxon, quant))))
}

#' @rdname el_drarefy
#' @export
el_rarefy <- function(x, depth, check = TRUE, sample = "sample", taxon = "taxon", quant = "quant") {
  x <- rarefy_input(x, depth, check, sample, taxon, quant)
  d <- as.data.frame(dplyr::summarise(x, k = dplyr::n(), .by = c(sample, total, quant)))
  c(tapply(d$k * presence_probability(d$total, d$quant, depth), d$sample, sum))
}

presence_probability <- function(total, quant, depth) {
  n <- pmin(depth, total)
  1 - ifelse(total - quant < n, 0, exp(lchoose(total - quant, n) - lchoose(total, n)))
}

rarefy_input <- function(x, depth, check, sample, taxon, quant) {
  if (!is.numeric(depth) || length(depth) != 1 || !is.finite(depth) || depth < 1 || depth != round(depth)) {
    stop("`depth` must be a single finite whole number, at least 1", call. = FALSE)
  }
  x <- duckplyr::as_duckdb_tibble(x, prudence = "stingy")
  x <- dplyr::select(x, sample = dplyr::all_of(sample), taxon = dplyr::all_of(taxon),
                     quant = dplyr::all_of(quant))
  if (check) {
    check_quantities(x, quant)
    fractional <- dplyr::summarise(x, n = sum(dplyr::if_else(round(quant) != quant, 1L, 0L), na.rm = TRUE))
    if (as.data.frame(fractional)$n > 0) stop("`", quant, "` must hold whole counts (integers)", call. = FALSE)
  }
  x <- dplyr::filter(x, quant > 0)
  x <- dplyr::mutate(x, total = sum(quant, na.rm = TRUE), .by = sample)
  if (check) {
    smallest <- as.data.frame(dplyr::summarise(x, smallest = min(total, na.rm = TRUE)))$smallest
    if (isTRUE(smallest < depth)) warning("some samples have a total below `depth` and are not rarefied", call. = FALSE)
  }
  x
}
