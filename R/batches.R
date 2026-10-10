utils::globalVariables("taxon")

# The query is written to a temporary Parquet file, then sorted and numbered by a
# DuckDB connection of our own, so that its memory limit and spill directory do
# not touch the user's duckplyr session. The kernel reads integer ids in batches
# that always hold whole taxa.
open_taxon_stream <- function(x, binary = FALSE, memory_limit = "250MB", threads = 1L) {
  x <- duckplyr::as_duckdb_tibble(x, prudence = "stingy")
  path <- tempfile(fileext = ".parquet")
  duckplyr::compute_parquet(x, path)

  con <- DBI::dbConnect(duckdb::duckdb())
  spill <- tempfile("ecolongr-spill")
  dir.create(spill)
  DBI::dbExecute(con, paste("SET memory_limit =", DBI::dbQuoteString(con, memory_limit)))
  DBI::dbExecute(con, sprintf("SET threads = %d", min(as.integer(threads), 4L)))
  DBI::dbExecute(con, paste("SET temp_directory =", DBI::dbQuoteString(con, spill)))
  DBI::dbExecute(con, "SET preserve_insertion_order = false")
  DBI::dbExecute(con, paste(
    "CREATE VIEW v AS SELECT sample, taxon, quant FROM read_parquet(", DBI::dbQuoteString(con, path), ")"
  ))

  structure(
    list(con = con, files = c(path, spill), sorted = tempfile(fileext = ".parquet"),
         value = if (binary) "1.0" else "quant"),
    class = "taxon_stream"
  )
}

close_taxon_stream <- function(stream) {
  DBI::dbDisconnect(stream$con, shutdown = TRUE)
  unlink(c(stream$files, stream$sorted), recursive = TRUE)
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
      SELECT CAST(sample AS VARCHAR) AS sample,
             coalesce(sum(%1$s) FILTER (WHERE quant <> 0), 0) AS total,
             coalesce(sum(%1$s * %1$s) FILTER (WHERE quant <> 0), 0) AS square
      FROM v GROUP BY sample", stream$value))
  )
}

stream_batches <- function(stream, samples, batch_rows, fun) {
  con <- stream$con
  duckdb::duckdb_register(con, "sample_ids", data.frame(
    sample = as.character(samples), id = seq_along(samples) - 1L
  ))
  on.exit(duckdb::duckdb_unregister(con, "sample_ids"), add = TRUE)

  sort_to_file <- function() {
    unlink(stream$sorted)
    DBI::dbExecute(con, paste("COPY (SELECT * FROM (",
      sprintf("SELECT s.id AS sample, CAST(dense_rank() OVER (ORDER BY v.taxon) AS INTEGER) AS taxon, %s AS quant", stream$value),
      "FROM v JOIN sample_ids s ON CAST(v.sample AS VARCHAR) = s.sample WHERE v.quant <> 0",
      ") ORDER BY taxon) TO", DBI::dbQuoteString(con, stream$sorted), "(FORMAT PARQUET)"))
  }
  # Each DuckDB thread needs tens of MB; with a small limit, one thread still works.
  tryCatch(sort_to_file(), error = function(e) {
    if (!grepl("Out of Memory", conditionMessage(e), fixed = TRUE)) stop(e)
    DBI::dbExecute(con, "SET threads = 1")
    sort_to_file()
  })

  res <- DBI::dbSendQuery(con, paste(
    "SELECT sample, taxon, quant FROM read_parquet(", DBI::dbQuoteString(con, stream$sorted), ")"))
  on.exit(DBI::dbClearResult(res), add = TRUE)

  emit <- function(v, from, to) {
    i <- from:to
    fun(data.frame(sample = v$sample[i], taxon = v$taxon[i], quant = v$quant[i]))
  }

  carry <- NULL
  repeat {
    rows <- DBI::dbFetch(res, n = batch_rows)
    done <- DBI::dbHasCompleted(res)
    v <- if (is.null(carry)) as.list(rows) else
      list(sample = c(carry$sample, rows$sample), taxon = c(carry$taxon, rows$taxon),
           quant = c(carry$quant, rows$quant))
    n <- length(v$taxon)
    if (n == 0) break
    if (done && n <= batch_rows) {
      emit(v, 1L, n)
      break
    }
    first_of_last <- sum(v$taxon < v$taxon[n]) + 1L
    if (first_of_last > 1L) emit(v, 1L, first_of_last - 1L)
    carry <- lapply(v, function(col) col[first_of_last:n])
    if (done) {
      emit(carry, 1L, length(carry$taxon))
      break
    }
  }
  invisible()
}

for_each_taxon_batch <- function(x, batch_rows, fun, samples = NULL) {
  stream <- open_taxon_stream(x)
  on.exit(close_taxon_stream(stream))
  if (is.null(samples)) samples <- sort(stream_summary(stream)$totals$sample)
  stream_batches(stream, samples, batch_rows, fun)
}

stream_taxa <- function(stream) {
  DBI::dbGetQuery(stream$con, "
    SELECT CAST(taxon AS VARCHAR) AS name, CAST(dense_rank() OVER (ORDER BY taxon) AS INTEGER) AS id
    FROM (SELECT DISTINCT taxon FROM v WHERE quant <> 0) ORDER BY id")$name
}
