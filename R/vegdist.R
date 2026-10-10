utils::globalVariables(c("sample_a", "sample_b", "quant_a", "quant_b"))

#' Bray-Curtis dissimilarity from a long table
#'
#' Computes the distance in the database from a table without zeros.
#' Only the sample pairs sharing at least one taxon are collected into R.
#'
#' @param x A data frame or lazy `duckplyr` table with one row per non-zero
#'   sample and taxon.
#' @param method Only `"bray"` is implemented.
#' @param samples Optional vector of all sample IDs. Samples in `samples` that
#'   are absent from `x` are treated as all zero.
#'   Without it, the samples that occur in `x` are used.
#' @param sample,taxon,quant Names of the columns in `x`.
#' @return A `dist` object.
#' @export
el_vegdist <- function(x, method = "bray", samples = NULL,
                       sample = "sample", taxon = "taxon", quant = "quant") {
  method <- match.arg(method, "bray")

  x <- duckplyr::as_duckdb_tibble(x, prudence = "stingy")
  x <- dplyr::select(x, sample = dplyr::all_of(sample), taxon = dplyr::all_of(taxon),
                     quant = dplyr::all_of(quant))

  totals <- as.data.frame(dplyr::summarise(x, total = sum(quant), .by = sample))
  totals$sample <- as.character(totals$sample)

  if (is.null(samples)) {
    samples <- sort(totals$sample)
  } else {
    samples <- as.character(samples)
    missing <- setdiff(totals$sample, samples)
    if (length(missing) > 0) {
      stop("Samples in `x` but not in `samples`: ", paste(utils::head(missing, 5), collapse = ", "))
    }
  }

  a <- dplyr::rename(x, sample_a = sample, quant_a = quant)
  b <- dplyr::rename(x, sample_b = sample, quant_b = quant)
  shared <- dplyr::inner_join(a, b, by = "taxon", relationship = "many-to-many")
  shared <- dplyr::filter(shared, sample_a < sample_b)
  shared <- dplyr::summarise(
    shared,
    shared_min = sum(dplyr::if_else(quant_a < quant_b, quant_a, quant_b)),
    .by = c(sample_a, sample_b)
  )
  shared <- as.data.frame(shared)

  total <- stats::setNames(totals$total, totals$sample)
  d <- matrix(1, length(samples), length(samples), dimnames = list(samples, samples))
  diag(d) <- 0

  i <- match(as.character(shared$sample_a), samples)
  j <- match(as.character(shared$sample_b), samples)
  value <- 1 - 2 * shared$shared_min / (total[as.character(shared$sample_a)] + total[as.character(shared$sample_b)])
  d[cbind(i, j)] <- value
  d[cbind(j, i)] <- value

  empty <- samples[!samples %in% totals$sample | samples %in% totals$sample[totals$total == 0]]
  d[empty, empty] <- NaN
  diag(d) <- 0

  out <- stats::as.dist(d)
  attr(out, "method") <- method
  out
}
