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
#' With no `depth`, the smallest sample total is used, so no sample is dropped.
#' A sample with a total below `depth` is dropped from the result, with a
#' warning that names it, and it is an error if no sample is left.
#' `vegan` instead treats such a sample as if `depth` were its total.
#' A sample whose total equals `depth` is kept.
#'
#' The probability depends only on the sample total and the count, so it is
#' computed once for each distinct pair of the two, in R with [lchoose()] as
#' `vegan` does. Only these pairs, not the table, are read into R.
#'
#' @inheritParams el_rrarefy
#' @param depth A single finite whole number, at least 1. The default, `NULL`,
#'   is the smallest total among the samples with a non-zero count, which takes
#'   one small query.
#' @param check If `TRUE`, run queries that stop on negative, missing, non-finite
#'   or non-whole counts, warn about samples below `depth` and stop if none is left.
#'   Use `FALSE` to stay lazy: small samples are still dropped, silently.
#' @return `el_drarefy()`: a lazy `duckplyr` table with the columns `sample`,
#'   `taxon` and `quant`, named as in the arguments. `el_rarefy()`: a numeric
#'   vector with one expected richness per sample, named by sample.
#' @example man-roxygen/lmo-counts.R
#' @examples
#' # Expected richness of each sample in a subsample of 5,000 reads.
#' expected <- el_rarefy(counts, depth = 5000)
#' observed <- el_specnumber(counts)
#' head(cbind(observed, expected = expected[names(observed)]))
#'
#' # The probability that each ASV is in such a subsample, as a lazy table.
#' counts |>
#'   el_drarefy(depth = 5000) |>
#'   head() |>
#'   as.data.frame()
#' @export
el_drarefy <- function(x, depth = NULL, check = TRUE, sample = "sample", taxon = "taxon", quant = "quant") {
  prepared <- rarefy_input(x, depth, check, sample, taxon, quant)
  x <- prepared$x
  depth <- prepared$depth
  pairs <- as.data.frame(dplyr::distinct(x, total, quant))
  pairs$p <- presence_probability(pairs$total, pairs$quant, depth)
  pairs <- duckplyr::as_duckdb_tibble(pairs, prudence = "stingy")
  out <- dplyr::inner_join(x, pairs, by = c("total", "quant"))
  out <- dplyr::select(out, sample, taxon, quant = p)
  dplyr::rename(out, dplyr::all_of(stats::setNames(c("sample", "taxon", "quant"), c(sample, taxon, quant))))
}

#' @rdname el_drarefy
#' @export
el_rarefy <- function(x, depth = NULL, check = TRUE, sample = "sample", taxon = "taxon", quant = "quant") {
  prepared <- rarefy_input(x, depth, check, sample, taxon, quant)
  x <- prepared$x
  depth <- prepared$depth
  d <- as.data.frame(dplyr::summarise(x, k = dplyr::n(), .by = c(sample, total, quant)))
  c(tapply(d$k * presence_probability(d$total, d$quant, depth), d$sample, sum))
}

presence_probability <- function(total, quant, depth) {
  1 - ifelse(total - quant < depth, 0, exp(lchoose(total - quant, depth) - lchoose(total, depth)))
}

rarefy_input <- function(x, depth, check, sample, taxon, quant) {
  if (!is.null(depth)) check_depth(depth)
  x <- duckplyr::as_duckdb_tibble(x, prudence = "stingy")
  x <- dplyr::select(x, sample = dplyr::all_of(sample), taxon = dplyr::all_of(taxon),
                     quant = dplyr::all_of(quant))
  if (check) check_counts(x, quant)
  x <- dplyr::filter(x, quant > 0)
  x <- dplyr::mutate(x, total = sum(quant, na.rm = TRUE), .by = sample)
  if (is.null(depth)) {
    depth <- as.data.frame(dplyr::summarise(x, smallest = min(total, na.rm = TRUE)))$smallest
    if (!is.finite(depth)) stop("`", quant, "` has no non-zero counts", call. = FALSE)
  }
  if (check) {
    small <- as.data.frame(dplyr::distinct(dplyr::filter(x, total < depth), sample))$sample
    kept <- as.data.frame(dplyr::summarise(dplyr::filter(x, total >= depth), n = dplyr::n()))$n
    if (kept == 0) stop("no sample has a total of at least `depth` (", format(depth), ")", call. = FALSE)
    if (length(small) > 0) warn_dropped(sort(as.character(small)))
  }
  list(x = dplyr::filter(x, total >= depth), depth = depth)
}

check_depth <- function(depth) {
  if (!is.numeric(depth) || length(depth) != 1 || !is.finite(depth) || depth < 1 || depth != round(depth)) {
    stop("`depth` must be a single finite whole number, at least 1", call. = FALSE)
  }
}

check_counts <- function(x, quant) {
  check_quantities(x, quant)
  fractional <- dplyr::summarise(x, n = sum(dplyr::if_else(round(quant) != quant, 1L, 0L), na.rm = TRUE))
  if (as.data.frame(fractional)$n > 0) stop("`", quant, "` must hold whole counts (integers)", call. = FALSE)
  invisible()
}
