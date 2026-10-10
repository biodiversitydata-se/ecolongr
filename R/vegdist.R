utils::globalVariables(c("sample_a", "sample_b", "quant_a", "quant_b"))

#' Dissimilarity from a long table
#'
#' Computes the distance in the database from a table without zeros.
#' Only the sample pairs sharing at least one taxon are collected into R.
#'
#' @param x A data frame or lazy `duckplyr` table with one row per non-zero
#'   sample and taxon. Rows with `quant == 0` are ignored.
#' @param method One of `"bray"`, `"jaccard"` and `"euclidean"`, as in
#'   [vegan::vegdist()].
#' @param binary If `TRUE`, use presence and absence instead of quantities.
#' @param samples Optional vector of all sample IDs, which also sets the order
#'   of the result. Samples in `samples` that are absent from `x` are treated as
#'   all zero. Without it, the samples that occur in `x` are used, in
#'   alphabetical order.
#' @param sample,taxon,quant Names of the columns in `x`.
#' @return A `dist` object.
#' @export
el_vegdist <- function(x, method = c("bray", "jaccard", "euclidean"), binary = FALSE,
                       samples = NULL, sample = "sample", taxon = "taxon", quant = "quant") {
  method <- match.arg(method)

  x <- duckplyr::as_duckdb_tibble(x, prudence = "stingy")
  x <- dplyr::select(x, sample = dplyr::all_of(sample), taxon = dplyr::all_of(taxon),
                     quant = dplyr::all_of(quant))

  checks <- as.data.frame(dplyr::summarise(
    x,
    missing = sum(dplyr::if_else(is.na(quant), 1L, 0L)),
    negative = sum(dplyr::if_else(quant < 0, 1L, 0L), na.rm = TRUE)
  ))
  if (isTRUE(checks$missing > 0)) stop("`", quant, "` has missing values")
  if (isTRUE(checks$negative > 0)) stop("`", quant, "` has negative values")

  x <- dplyr::filter(x, quant > 0)
  if (binary) x <- dplyr::mutate(x, quant = 1)

  totals <- as.data.frame(dplyr::summarise(
    x, total = sum(quant), square = sum(quant * quant), .by = sample
  ))
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

  n <- length(samples)
  total <- square <- numeric(n)
  at <- match(totals$sample, samples)
  total[at] <- totals$total
  square[at] <- totals$square

  a <- dplyr::rename(x, sample_a = sample, quant_a = quant)
  b <- dplyr::rename(x, sample_b = sample, quant_b = quant)
  shared <- dplyr::filter(
    dplyr::inner_join(a, b, by = "taxon", relationship = "many-to-many"),
    sample_a < sample_b
  )
  shared <- if (method == "euclidean") {
    dplyr::summarise(shared, overlap = sum(quant_a * quant_b), .by = c(sample_a, sample_b))
  } else {
    dplyr::summarise(
      shared,
      overlap = sum(dplyr::if_else(quant_a < quant_b, quant_a, quant_b)),
      .by = c(sample_a, sample_b)
    )
  }
  shared <- as.data.frame(shared)
  i <- match(as.character(shared$sample_a), samples)
  j <- match(as.character(shared$sample_b), samples)

  if (method == "euclidean") {
    d <- sqrt(outer(square, square, "+"))
    value <- sqrt(pmax(square[i] + square[j] - 2 * shared$overlap, 0))
  } else {
    d <- matrix(1, n, n)
    value <- 1 - 2 * shared$overlap / (total[i] + total[j])
    if (method == "jaccard") value <- 2 * value / (1 + value)
    d[outer(total == 0, total == 0, "&")] <- NaN
  }
  d[cbind(i, j)] <- value
  d[cbind(j, i)] <- value
  diag(d) <- 0
  dimnames(d) <- list(samples, samples)

  out <- stats::as.dist(d)
  attr(out, "method") <- method
  out
}
