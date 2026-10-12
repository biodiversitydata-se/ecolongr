#' Dissimilarity from a long table
#'
#' Computes the distance from a table without zeros.
#' The query is run by DuckDB and written to a temporary Parquet file sorted by
#' taxon. A compiled kernel then reads it in batches and adds up the shared
#' quantities for every pair of samples, so the memory used is the result
#' (one number per pair of samples) plus one batch.
#'
#' @param x A data frame or lazy `duckplyr` table with one row per non-zero
#'   sample and taxon. Rows with `quant == 0` are ignored, but a sample that
#'   has only such rows is kept as an all-zero sample. Negative values are
#'   allowed only for `method = "euclidean"` without `binary`, for example after
#'   `el_decostand(method = "rclr")`.
#' @param method One of `"bray"`, `"jaccard"`, `"euclidean"`, `"aitchison"`,
#'   `"manhattan"` and `"canberra"`, as in [vegan::vegdist()]. `"aitchison"` is
#'   the Euclidean distance between clr-transformed samples, computed from the
#'   counts. `"canberra"` is `NA`, as in `vegan`, between two all-zero samples.
#'   `"robust.aitchison"` is the Euclidean distance on the robust clr with the
#'   zeros filled by matrix completion (rank `ropt`, as in `vegan`). It needs
#'   counts, and no all-zero samples.
#' @param binary If `TRUE`, use presence and absence instead of quantities.
#' @param pseudocount For `method = "aitchison"` only: added to the counts before
#'   the clr. A positive number, `"inverse_depth"` for `1 / (the sample's total
#'   count)`, or a positive numeric vector named by sample. The unit is that of
#'   the counts, so give raw counts, not relative abundances.
#'   `"inverse_depth"` uses the totals of the table passed in, so after filtering
#'   taxa it is not the depth of the whole library: compute the depths first and
#'   pass them as a named vector.
#' @param ropt For `method = "robust.aitchison"` only: the rank of the matrix
#'   completion, 3 as in `vegan`. A whole number below the number of samples and
#'   of taxa.
#' @param samples Optional vector of all sample IDs, which also sets the order
#'   of the result. Samples in `samples` that are absent from `x` are treated as
#'   all zero. Without it, the samples that occur in `x` are used, in
#'   alphabetical order.
#' @details
#' Euclidean distance on `el_decostand(method = "rclr")` output reads absent
#' cells as 0. That is not `vegan`'s `"robust.aitchison"` when the table has
#' zeros, because `vegan` first fills them by matrix completion; see
#' [el_decostand()].
#'
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
el_vegdist <- function(x, method = c("bray", "jaccard", "euclidean", "aitchison", "manhattan", "canberra", "robust.aitchison"), binary = FALSE,
                       pseudocount = 1, ropt = 3, samples = NULL, batch_rows = 1e6, threads = NULL, memory_limit = NULL,
                       sample = "sample", taxon = "taxon", quant = "quant") {
  method <- match.arg(method)
  if (method %in% c("aitchison", "robust.aitchison") && binary) {
    stop("`binary` is not available for method = \"", method, "\"")
  }
  if (method == "aitchison") check_pseudocount_form(pseudocount)
  if (method != "robust.aitchison" && !missing(ropt)) {
    stop("`ropt` is only used with `method = \"robust.aitchison\"`")
  }
  if (!is.numeric(batch_rows) || length(batch_rows) != 1 || is.na(batch_rows) || batch_rows < 1) {
    stop("`batch_rows` must be a positive number")
  }

  threads <- resolve_threads(threads)
  memory_limit <- resolve_memory_limit(memory_limit)

  x <- duckplyr::as_duckdb_tibble(x, prudence = "stingy")
  x <- dplyr::select(x, sample = dplyr::all_of(sample), taxon = dplyr::all_of(taxon),
                     quant = dplyr::all_of(quant))

  stream <- open_taxon_stream(x, binary, memory_limit = memory_limit, threads = threads)
  on.exit(close_taxon_stream(stream))

  summary <- stream_summary(stream)
  if (summary$missing > 0) stop("`", quant, "` has missing values")
  if (summary$negative > 0 && (method != "euclidean" || binary)) stop("`", quant, "` has negative values")
  if (method == "robust.aitchison") {
    fit <- rclr_completion(stream, summary, if (is.null(samples)) NULL else as.character(samples), threads, ropt)
    out <- stats::dist(fit$coordinates)
    attr(out, "Labels") <- fit$samples
    attr(out, "method") <- method
    return(out)
  }
  if (method == "aitchison") {
    stream_clr(stream, resolve_pseudocount(pseudocount, summary$totals))
    summary <- stream_summary(stream)
  }

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
  total <- square <- cells <- numeric(n)
  at <- match(totals$sample, samples)
  total[at] <- totals$total
  square[at] <- totals$square
  cells[at] <- totals$cells

  acc <- numeric(n * (n - 1) / 2)
  shared <- numeric(n * (n - 1) / 2)
  product <- method %in% c("euclidean", "aitchison")
  canberra <- method == "canberra"
  stream_batches(stream, samples, batch_rows, function(rows) {
    accumulate_pairs(acc, rows$sample, rows$taxon, rows$quant, n, product, threads = threads,
                     canberra = canberra, count = if (canberra) shared else numeric(0))
  })
  finalize_dist(acc, total, square, n,
                match(method, c("bray", "jaccard", "euclidean", "aitchison", "manhattan", "canberra")) - 1L,
                threads = threads, ntaxa = if (method == "aitchison") length(stream_taxa(stream)) else 1,
                nnz = cells, shared = if (canberra) shared else numeric(0))

  structure(acc, Size = n, Labels = samples, Diag = FALSE, Upper = FALSE,
            method = method, class = "dist")
}
