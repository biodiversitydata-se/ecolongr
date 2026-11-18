#' Rarefy samples to a common depth
#'
#' Counterpart of [vegan::rrarefy()]: a random subsample, without replacement,
#' of `depth` individuals from each sample. Samples with a total at or below
#' `depth` are returned unchanged, with a warning if the total is below.
#'
#' The result is random, so call [set.seed()] first for a reproducible one.
#' It is not the same random sample as `vegan`'s for the same seed: `vegan`
#' draws the individuals one at a time, and `el_rrarefy()` draws the count of
#' each taxon from the hypergeometric distribution of what is left. Both draw
#' from the same distribution of counts (the multivariate hypergeometric).
#'
#' The table is read into R and copied several times, so this is not a lazy
#' step: allow a few times the memory of the input. Counts or totals of about
#' 2^31 and above make each draw slow.
#'
#' @param x A data frame or lazy `duckplyr` table with one row per non-zero
#'   sample and taxon, holding whole counts. Rows with `quant == 0` are dropped.
#' @param depth A single finite whole number, at least 1.
#' @param sample,taxon,quant Names of the columns in `x`. The result has only
#'   these three columns, with the same names.
#' @return A `duckplyr` table holding the result in memory, with the samples in
#'   the order they first appear in `x` and without the taxa that were not drawn.
#' @export
el_rrarefy <- function(x, depth, sample = "sample", taxon = "taxon", quant = "quant") {
  if (!is.numeric(depth) || length(depth) != 1 || !is.finite(depth) || depth < 1 || depth != round(depth)) {
    stop("`depth` must be a single finite whole number, at least 1", call. = FALSE)
  }
  x <- duckplyr::as_duckdb_tibble(x, prudence = "stingy")
  d <- as.data.frame(dplyr::select(x, sample = dplyr::all_of(sample), taxon = dplyr::all_of(taxon),
                                   quant = dplyr::all_of(quant)))
  d$quant <- as.numeric(d$quant)
  if (anyNA(d$quant)) stop("`", quant, "` has missing values", call. = FALSE)
  if (any(d$quant < 0)) stop("`", quant, "` has negative values", call. = FALSE)
  if (any(!is.finite(d$quant) | d$quant != round(d$quant))) {
    stop("`", quant, "` must hold finite whole counts (integers)", call. = FALSE)
  }
  d <- d[d$quant > 0, ]
  id <- match(d$sample, unique(d$sample))
  o <- order(id)
  d <- d[o, ]
  id <- id[o]
  if (any(rowsum(d$quant, id) < depth)) warning("some samples have a total below `depth` and are not rarefied", call. = FALSE)
  d$quant <- rarefy_runs(d$quant, c(0L, which(diff(id) != 0L), nrow(d)), depth)
  d <- d[d$quant > 0, ]
  rownames(d) <- NULL
  names(d) <- c(sample, taxon, quant)
  duckplyr::as_duckdb_tibble(d, prudence = "stingy")
}
