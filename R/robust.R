# The robust clr with matrix completion, as vegan::decostand(x, "rclr") with impute = TRUE:
# the clr of the non-zero cells is completed to a rank-ropt matrix with OptSpace, and the
# result is centred on rows and columns. The completed matrix is X S Y' with small factors,
# so it is never made dense: after centring, M = Xc S Yc', and with Yc = Qb Rb the rows of
# Xc S Rb' have the same distances as the rows of M.
rclr_completion <- function(stream, summary, samples, threads, ropt = 3) {
  if (any(summary$totals$cells == 0) ||
      (!is.null(samples) && !setequal(samples, as.character(summary$totals$sample)))) {
    stop("all-zero samples are not supported with the completed rclr: they have no observed cells to complete from", call. = FALSE)
  }
  if (is.null(samples)) samples <- sort(as.character(summary$totals$sample))
  stream <- stream_rclr(stream)
  taxa <- stream_taxa(stream)
  if (!is.numeric(ropt) || length(ropt) != 1 || is.na(ropt) || ropt < 1 || ropt != round(ropt) ||
      ropt >= min(length(samples), length(taxa))) {
    stop("`ropt` must be a positive whole number below the number of samples and of taxa", call. = FALSE)
  }

  rows <- list()
  stream_batches(stream, samples, 1e6, function(b) rows[[length(rows) + 1L]] <<- b)
  i <- unlist(lapply(rows, `[[`, "sample"), use.names = FALSE)
  j <- unlist(lapply(rows, `[[`, "taxon"), use.names = FALSE) - 1L
  v <- unlist(lapply(rows, `[[`, "quant"), use.names = FALSE)
  rm(rows)

  fit <- optspace_triplets(i, j, v, length(samples), length(taxa), r = ropt, threads = threads)
  Xc <- sweep(fit$X, 2, colMeans(fit$X))
  Yc <- sweep(fit$Y, 2, colMeans(fit$Y))
  qa <- qr(Xc)
  qb <- qr(Yc)
  Ra <- qr.R(qa)
  Rb <- qr.R(qb)
  list(
    samples = samples,
    taxa = taxa,
    coordinates = Xc %*% fit$S %*% t(Rb),
    Qa = qr.Q(qa),
    Qb = qr.Q(qb),
    core = Ra %*% fit$S %*% t(Rb)
  )
}
