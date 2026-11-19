#' Rarefy samples to a common depth
#'
#' Counterpart of [vegan::rrarefy()]: a random subsample, without replacement,
#' of `depth` individuals from each sample. With no `depth`, the smallest
#' sample total is used, so no sample is dropped.
#'
#' A sample with a total below `depth` cannot be rarefied and is dropped from the
#' result, with a warning that names it; it is an error if no sample is left.
#' `vegan` returns such a sample unchanged, which leaves samples of different
#' depths in the result. A sample whose total equals `depth` is kept as it is.
#'
#' The result is random, so call [set.seed()] first for a reproducible one.
#' The rows are sorted by sample, taxon and count before anything is drawn, so
#' the same seed gives the same table however the rows of `x` are ordered, which
#' for a lazy table is not under your control.
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
#' @param depth A single finite whole number, at least 1. The default, `NULL`,
#'   is the smallest total among the samples with a non-zero count.
#' @param sample,taxon,quant Names of the columns in `x`. The result has only
#'   these three columns, with the same names.
#' @return A `duckplyr` table holding the result in memory, sorted by sample and
#'   taxon and without the taxa that were not drawn.
#' @example man-roxygen/lmo-counts.R
#' @examples
#' set.seed(1)
#' subsample <- el_rrarefy(counts, depth = 5000)
#'
#' # Every sample now has exactly 5,000 reads, but different samples keep different ASVs.
#' subsample |> summarise(reads = sum(quant), .by = sample) |> head() |> as.data.frame()
#' el_specnumber(subsample) |> head()
#' @export
el_rrarefy <- function(x, depth = NULL, sample = "sample", taxon = "taxon", quant = "quant") {
  if (!is.null(depth)) check_depth(depth)
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
  d <- d[order(d$sample, d$taxon, d$quant, method = "radix"), ]
  samples <- unique(d$sample)
  total <- c(rowsum(d$quant, match(d$sample, samples)))
  if (length(total) == 0 && is.null(depth)) stop("`", quant, "` has no non-zero counts", call. = FALSE)
  if (is.null(depth)) depth <- min(total)
  if (all(total < depth)) stop("no sample has a total of at least `depth` (", format(depth), ")", call. = FALSE)
  if (any(total < depth)) {
    warn_dropped(samples[total < depth])
    d <- d[d$sample %in% samples[total >= depth], ]
  }
  id <- match(d$sample, unique(d$sample))
  d$quant <- rarefy_runs(d$quant, c(0L, which(diff(id) != 0L), nrow(d)), depth)
  d <- d[d$quant > 0, ]
  rownames(d) <- NULL
  names(d) <- c(sample, taxon, quant)
  duckplyr::as_duckdb_tibble(d, prudence = "stingy")
}

warn_dropped <- function(dropped) {
  shown <- paste(utils::head(dropped, 5), collapse = ", ")
  if (length(dropped) > 5) shown <- paste0(shown, ", ...")
  warning(length(dropped), if (length(dropped) == 1) " sample" else " samples",
          " with a total below `depth` ", if (length(dropped) == 1) "was" else "were", " dropped: ", shown,
          call. = FALSE)
}
