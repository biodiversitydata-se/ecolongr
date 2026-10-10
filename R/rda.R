#' Principal components of a long table
#'
#' Unconstrained [vegan::rda()] on a table without zeros. The matrix is
#' centred implicitly, so it is never made dense, and only the first `k` axes
#' are returned. Compose it with [el_decostand()] for a transformation.
#'
#' @param x A data frame or lazy `duckplyr` table with one row per non-zero
#'   sample and taxon. Negative values are allowed.
#' @param k Number of axes to return, at most the rank, `min(n - 1, taxa)`.
#' @param scale If `TRUE`, divide every taxon by its standard deviation, as in
#'   `vegan::rda(scale = TRUE)`.
#' @param samples Optional vector of all sample IDs. Samples in `samples` that
#'   are absent from `x` are all-zero rows and count in the means and the total
#'   inertia.
#' @param transform `"none"` uses `x` as it is. `"clr"` takes raw counts and
#'   centred log-ratio transforms them, with every zero cell read as
#'   `log(pseudocount)`, without making the matrix dense. `"rclr"` takes raw
#'   counts, fills the zeros of the robust clr by rank-3 matrix completion as in
#'   `vegan::decostand(x, "rclr")`, and returns at most three axes.
#' @param pseudocount For `transform = "clr"` only: a positive number, `"inverse_depth"`
#'   for `1 / (the sample's total count)`, or a positive numeric vector named by
#'   sample. The unit is that of the counts.
#' @param solver `"exact"` uses a dense SVD and is meant for small tables,
#'   `"irlba"` a sparse truncated SVD, and `"auto"` chooses by size.
#' @param sample,taxon,quant Names of the columns in `x`.
#' @return An object of class `rda` and `cca` with the first `k` axes. Use
#'   [vegan::scores()] and [vegan::eigenvals()] on it. Functions that need the
#'   full fitted matrix, such as `predict()` and `anova()`, do not work.
#' @export
el_rda <- function(x, k = 5, scale = FALSE, samples = NULL, solver = c("auto", "exact", "irlba"),
                   transform = c("none", "clr", "rclr"), pseudocount = 1,
                   sample = "sample", taxon = "taxon", quant = "quant") {
  call <- match.call()
  solver <- match.arg(solver)
  transform <- match.arg(transform)
  if (transform != "clr" && !missing(pseudocount)) {
    stop("`pseudocount` is only used with `transform = \"clr\"`")
  }
  if (transform == "rclr" && scale) stop("`scale = TRUE` is not available with `transform = \"rclr\"`")
  if (transform == "clr") check_pseudocount_form(pseudocount)
  if (!is.numeric(k) || length(k) != 1 || is.na(k) || k < 1 || k != round(k)) {
    stop("`k` must be a positive whole number")
  }

  x <- duckplyr::as_duckdb_tibble(x, prudence = "stingy")
  x <- dplyr::select(x, sample = dplyr::all_of(sample), taxon = dplyr::all_of(taxon),
                     quant = dplyr::all_of(quant))

  stream <- open_taxon_stream(x, memory_limit = resolve_memory_limit(NULL), threads = resolve_threads(NULL))
  on.exit(close_taxon_stream(stream))

  summary <- stream_summary(stream)
  if (summary$missing > 0) stop("`", quant, "` has missing values")
  if (transform == "rclr") {
    if (summary$negative > 0) stop("`", quant, "` has negative values")
    fit <- rclr_completion(stream, summary, if (is.null(samples)) NULL else as.character(samples), resolve_threads(NULL))
    n <- length(fit$samples)
    k <- min(k, ncol(fit$core))
    sv <- svd(fit$core)
    u <- fit$Qa %*% sv$u[, seq_len(k), drop = FALSE]
    vv <- fit$Qb %*% sv$v[, seq_len(k), drop = FALSE]
    axes <- paste0("PC", seq_len(k))
    dimnames(u) <- list(fit$samples, axes)
    dimnames(vv) <- list(fit$taxa, axes)
    tot_chi <- sum(sv$d^2) / (n - 1)
    colsum <- stats::setNames(sqrt(rowSums((fit$Qb %*% t(fit$core))^2) / (n - 1)), fit$taxa)
    return(structure(
      list(colsum = colsum, tot.chi = tot_chi, Ybar = NULL, method = "rda", call = call,
           pCCA = NULL, CCA = NULL,
           CA = list(eig = stats::setNames(sv$d[seq_len(k)]^2 / (n - 1), axes), poseig = NULL, u = u, v = vv,
                     rank = k, tot.chi = tot_chi),
           inertia = "variance",
           ecolongr = list(solver = "optspace", n = n, taxa = length(fit$taxa), k = k)),
      class = c("rda", "cca")
    ))
  }
  if (transform == "clr") {
    if (summary$negative > 0) stop("`", quant, "` has negative values")
    stream_clr(stream, resolve_pseudocount(pseudocount, summary$totals))
    summary <- stream_summary(stream)
  }
  sample_names <- as.character(summary$totals$sample)
  if (is.null(samples)) {
    samples <- sort(sample_names)
  } else {
    samples <- as.character(samples)
    unlisted <- setdiff(sample_names, samples)
    if (length(unlisted) > 0) {
      stop("Samples in `x` but not in `samples`: ", paste(utils::head(unlisted, 5), collapse = ", "))
    }
  }
  taxa <- stream_taxa(stream)

  rows <- list()
  stream_batches(stream, samples, 1e6, function(b) rows[[length(rows) + 1L]] <<- b)
  i <- unlist(lapply(rows, `[[`, "sample"), use.names = FALSE)
  j <- unlist(lapply(rows, `[[`, "taxon"), use.names = FALSE)
  v <- unlist(lapply(rows, `[[`, "quant"), use.names = FALSE)
  n <- length(samples)
  p <- length(taxa)
  A <- Matrix::sparseMatrix(i = i + 1L, j = j, x = v, dims = c(n, p))

  clr <- transform == "clr"
  if (clr) {
    row_sum <- numeric(n)
    row_sum[match(as.character(summary$totals$sample), samples)] <- summary$totals$total
    u <- row_sum / p
    centre <- Matrix::colSums(A) / n - mean(u)
    sum_squares <- Matrix::colSums(A^2) - 2 * Matrix::colSums(Matrix::Diagonal(x = u) %*% A) + sum(u^2)
  } else {
    centre <- Matrix::colSums(A) / n
    sum_squares <- Matrix::colSums(A^2)
  }
  variance <- pmax((sum_squares - n * centre^2) / (n - 1), 0)
  sd <- sqrt(variance)
  if (scale) {
    if (any(sd == 0)) stop("`scale = TRUE` needs taxa that are not constant; ", sum(sd == 0), " constant taxa found")
    tot_chi <- p
    colsum <- rep(1, p)
  } else {
    tot_chi <- sum(variance)
    colsum <- sd
  }

  rank <- min(n - 1L, p)
  k <- min(k, rank)
  if (solver == "auto") solver <- if (as.numeric(n) * p <= 2e5 || k >= min(n, p) - 1L) "exact" else "irlba"
  if (solver == "irlba" && k >= min(n, p)) solver <- "exact"
  usv <- if (solver == "exact") {
    z <- as.matrix(A)
    if (clr) z <- z - u
    z <- sweep(z, 2, colMeans(z))
    if (scale) z <- sweep(z, 2, sd, "/")
    s <- svd(z, nu = k, nv = k)
    list(d = s$d[seq_len(k)], u = s$u, v = s$v)
  } else if (clr) {
    op <- methods::new("clr_operator", A, u = u, cvec = centre, scale = if (scale) sd else rep(1, p))
    irlba::irlba(op, nv = k, tol = 1e-10, maxit = 5000)
  } else {
    irlba::irlba(A, nv = k, center = centre, scale = if (scale) sd else NULL, tol = 1e-10, maxit = 5000)
  }

  flip <- sign(vapply(seq_len(k), function(a) usv$u[which.max(abs(usv$u[, a])), a], numeric(1)))
  u <- sweep(usv$u, 2, flip, "*")
  vv <- sweep(usv$v, 2, flip, "*")
  axes <- paste0("PC", seq_len(k))
  dimnames(u) <- list(samples, axes)
  dimnames(vv) <- list(taxa, axes)
  eig <- stats::setNames(usv$d^2 / (n - 1), axes)
  names(colsum) <- taxa

  structure(
    list(colsum = colsum, tot.chi = tot_chi, Ybar = NULL, method = "rda", call = call,
         pCCA = NULL, CCA = NULL,
         CA = list(eig = eig, poseig = NULL, u = u, v = vv, rank = k, tot.chi = tot_chi),
         inertia = "variance",
         ecolongr = list(solver = solver, n = n, taxa = p, k = k)),
    class = c("rda", "cca")
  )
}
