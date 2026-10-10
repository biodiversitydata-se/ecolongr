#' Dissimilarity from a long table
#'
#' Computes the distance from a table without zeros.
#' The query is run by DuckDB and written to a temporary Parquet file sorted by
#' taxon. A compiled kernel then reads it in batches and adds up the shared
#' quantities for every pair of samples, so the memory used is the result
#' (one number per pair of samples) plus one batch.
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
#' @param batch_rows Approximate number of table rows held in R at a time.
#'   A batch always holds whole taxa.
#' @param sample,taxon,quant Names of the columns in `x`.
#' @return A `dist` object.
#' @export
el_vegdist <- function(x, method = c("bray", "jaccard", "euclidean"), binary = FALSE,
                       samples = NULL, batch_rows = 5e6,
                       sample = "sample", taxon = "taxon", quant = "quant") {
  method <- match.arg(method)
  if (!is.numeric(batch_rows) || length(batch_rows) != 1 || is.na(batch_rows) || batch_rows < 1) {
    stop("`batch_rows` must be a positive number")
  }

  x <- duckplyr::as_duckdb_tibble(x, prudence = "stingy")
  x <- dplyr::select(x, sample = dplyr::all_of(sample), taxon = dplyr::all_of(taxon),
                     quant = dplyr::all_of(quant))

  stream <- open_taxon_stream(x, binary)
  on.exit(close_taxon_stream(stream))

  summary <- stream_summary(stream)
  if (summary$missing > 0) stop("`", quant, "` has missing values")
  if (summary$negative > 0) stop("`", quant, "` has negative values")

  totals <- summary$totals
  totals$sample <- as.character(totals$sample)
  if (is.null(samples)) {
    samples <- sort(totals$sample)
  } else {
    samples <- as.character(samples)
    unlisted <- setdiff(totals$sample, samples)
    if (length(unlisted) > 0) {
      stop("Samples in `x` but not in `samples`: ", paste(utils::head(unlisted, 5), collapse = ", "))
    }
  }

  n <- length(samples)
  total <- square <- numeric(n)
  at <- match(totals$sample, samples)
  total[at] <- totals$total
  square[at] <- totals$square

  acc <- numeric(n * (n - 1) / 2)
  product <- method == "euclidean"
  stream_batches(stream, batch_rows, function(rows) {
    taxon_id <- cumsum(c(TRUE, rows$taxon[-1] != rows$taxon[-nrow(rows)]))
    accumulate_pairs(acc, match(rows$sample, samples) - 1L, as.integer(taxon_id),
                     rows$quant, n, product)
  })
  finalize_dist(acc, total, square, n, match(method, c("bray", "jaccard", "euclidean")) - 1L)

  structure(acc, Size = n, Labels = samples, Diag = FALSE, Upper = FALSE,
            method = method, class = "dist")
}
