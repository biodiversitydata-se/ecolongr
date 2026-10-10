utils::globalVariables("taxon")

# The query is written once to a temporary Parquet file sorted by taxon and then
# read back through a connection of our own, so the rows a taxon spans are
# contiguous and only one batch is in R at a time.
open_taxon_stream <- function(x, binary = FALSE) {
  x <- duckplyr::as_duckdb_tibble(x, prudence = "stingy")
  path <- tempfile(fileext = ".parquet")
  duckplyr::compute_parquet(dplyr::arrange(x, taxon), path)
  con <- DBI::dbConnect(duckdb::duckdb())
  value <- if (binary) "1.0" else "quant"
  DBI::dbExecute(con, sprintf("CREATE VIEW v AS SELECT sample, taxon, quant FROM read_parquet('%s')", path))
  structure(list(con = con, path = path, value = value), class = "taxon_stream")
}

close_taxon_stream <- function(stream) {
  DBI::dbDisconnect(stream$con, shutdown = TRUE)
  unlink(stream$path)
  invisible()
}

stream_summary <- function(stream) {
  bad <- DBI::dbGetQuery(stream$con, "
    SELECT count(*) FILTER (WHERE quant IS NULL) AS missing,
           count(*) FILTER (WHERE quant < 0) AS negative
    FROM v")
  list(
    missing = bad$missing,
    negative = bad$negative,
    totals = DBI::dbGetQuery(stream$con, sprintf("
      SELECT sample, sum(%1$s) AS total, sum(%1$s * %1$s) AS square
      FROM v WHERE quant > 0 GROUP BY sample", stream$value))
  )
}

stream_batches <- function(stream, batch_rows, fun) {
  res <- DBI::dbSendQuery(stream$con, sprintf(
    "SELECT sample, taxon, %s AS quant FROM v WHERE quant > 0", stream$value))
  on.exit(DBI::dbClearResult(res), add = TRUE)

  carry <- NULL
  repeat {
    rows <- DBI::dbFetch(res, n = batch_rows)
    done <- DBI::dbHasCompleted(res)
    rows <- rbind(carry, rows)
    if (nrow(rows) == 0) break
    if (done && nrow(rows) <= batch_rows) {
      fun(rows)
      break
    }
    last <- rows$taxon == rows$taxon[nrow(rows)]
    if (!all(last)) fun(rows[!last, , drop = FALSE])
    carry <- rows[last, , drop = FALSE]
    if (done) {
      fun(carry)
      break
    }
  }
  invisible()
}

for_each_taxon_batch <- function(x, batch_rows, fun) {
  stream <- open_taxon_stream(x)
  on.exit(close_taxon_stream(stream))
  stream_batches(stream, batch_rows, fun)
}
