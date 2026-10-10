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
#' @param threads Number of threads for the pair kernel. Defaults to the option
#'   `ecolongr.threads`, or half of the cores. The result does not depend on it.
#' @param memory_limit Memory DuckDB may use for sorting the table, such as
#'   `"250MB"`; it spills to disk beyond that. Defaults to the option
#'   `ecolongr.memory_limit`, or `"250MB"`. This is separate from the memory
#'   used by R for the result and the current batch.
#' @param sample,taxon,quant Names of the columns in `x`.
#' @return A `dist` object.
#' @export
el_vegdist <- function(x, method = c("bray", "jaccard", "euclidean"), binary = FALSE,
                       samples = NULL, batch_rows = 1e6, threads = NULL, memory_limit = NULL,
                       sample = "sample", taxon = "taxon", quant = "quant") {
  method <- match.arg(method)
  if (!is.numeric(batch_rows) || length(batch_rows) != 1 || is.na(batch_rows) || batch_rows < 1) {
    stop("`batch_rows` must be a positive number")
  }

  if (is.null(threads)) threads <- getOption("ecolongr.threads", max(1L, parallel::detectCores() %/% 2L))
  if (!is.numeric(threads) || length(threads) != 1 || is.na(threads) || threads < 1 || threads != round(threads)) {
    stop("`threads` must be a positive whole number")
  }
  threads <- as.integer(threads)

  if (is.null(memory_limit)) memory_limit <- getOption("ecolongr.memory_limit", "250MB")
  if (!is.character(memory_limit) || length(memory_limit) != 1 ||
      !grepl("^[0-9]+(\\.[0-9]+)? ?[KMGT]i?B$", memory_limit)) {
    stop("`memory_limit` must be a single size such as \"250MB\"")
  }

  x <- duckplyr::as_duckdb_tibble(x, prudence = "stingy")
  x <- dplyr::select(x, sample = dplyr::all_of(sample), taxon = dplyr::all_of(taxon),
                     quant = dplyr::all_of(quant))

  stream <- open_taxon_stream(x, binary, memory_limit = memory_limit, threads = threads)
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
  stream_batches(stream, samples, batch_rows, function(rows) {
    accumulate_pairs(acc, rows$sample, rows$taxon, rows$quant, n, product, threads = threads)
  })
  finalize_dist(acc, total, square, n, match(method, c("bray", "jaccard", "euclidean")) - 1L, threads = threads)

  structure(acc, Size = n, Labels = samples, Diag = FALSE, Upper = FALSE,
            method = method, class = "dist")
}
