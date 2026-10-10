utils::globalVariables(c("sample", "quant", "minpos", "fractional"))

#' Standardise a long table within samples
#'
#' Lazy counterparts of the [vegan::decostand()] methods that need only the
#' non-zero cells. The result is a lazy `duckplyr` table that can go straight
#' into [el_vegdist()] or further `duckplyr` verbs.
#'
#' * `"total"`: quantities divided by the sample total.
#' * `"hellinger"`: square root of `"total"`.
#' * `"pa"`: 1 for every row.
#' * `"log"`: `log_b(x) + 1` as in `vegan`. If any value in the whole table is
#'   not a whole number, all values are first divided by the smallest positive
#'   value of the table, without a warning.
#' * `"rclr"`: `log(x)` minus the mean of `log(x)` over the sample's rows. This
#'   is `vegan`'s `"rclr"` with `impute = FALSE`; zeros stay absent where
#'   `vegan` returns `NA`. See the section on zeros below.
#'
#' @param x A data frame or lazy `duckplyr` table with one row per non-zero
#'   sample and taxon. Rows with `quant == 0` are dropped.
#' @param method One of the methods above.
#' @param logbase Base of the logarithm for `method = "log"`.
#' @param check If `TRUE`, run one query to stop on negative or missing
#'   quantities. Use `FALSE` to stay entirely lazy.
#' @param sample,taxon,quant Names of the columns in `x`. The result has only
#'   these three columns, with the same names.
#' @section Zeros and `rclr`:
#' `vegan`'s default, `decostand(x, "rclr")`, fills the missing cells by matrix
#' completion and then centres rows and columns, which is a global step that
#' changes the non-zero cells too. `el_decostand()` does not: zero cells stay
#' absent, and a later step that reads absent cells as 0 (such as
#' `el_vegdist(method = "euclidean")`) puts those taxa at the sample's
#' geometric mean. The result is therefore not `vegan`'s
#' `vegdist(method = "robust.aitchison")` when the table has zeros; the two
#' agree when it has none.
#'
#' @return A lazy `duckplyr` table.
#' @export
el_decostand <- function(x, method = c("total", "hellinger", "log", "pa", "rclr"),
                         logbase = 2, check = TRUE,
                         sample = "sample", taxon = "taxon", quant = "quant") {
  if (identical(method, "clr")) {
    stop("`clr` makes every cell non-zero, so it has no long-table form; use ",
         "el_vegdist(method = \"aitchison\") or el_rda(transform = \"clr\")", call. = FALSE)
  }
  method <- match.arg(method)
  if (method == "log" && (!is.numeric(logbase) || length(logbase) != 1 || is.na(logbase) || logbase <= 1)) {
    stop("`logbase` must be a number greater than 1")
  }

  x <- duckplyr::as_duckdb_tibble(x, prudence = "stingy")
  x <- dplyr::select(x, sample = dplyr::all_of(sample), taxon = dplyr::all_of(taxon),
                     quant = dplyr::all_of(quant))
  if (check) check_quantities(x, quant)
  x <- dplyr::filter(x, quant > 0)

  x <- switch(method,
    total = dplyr::mutate(x, quant = quant / sum(quant, na.rm = TRUE), .by = sample),
    hellinger = hellinger_transform(x),
    pa = dplyr::mutate(x, quant = 1),
    log = log_transform(x, logbase),
    rclr = dplyr::mutate(x, quant = log(quant) - mean(log(quant), na.rm = TRUE), .by = sample)
  )

  dplyr::rename(x, dplyr::all_of(stats::setNames(c("sample", "taxon", "quant"), c(sample, taxon, quant))))
}

# duckplyr cannot translate sqrt(), so the relative abundances go through a
# temporary Parquet file and the square root is taken in SQL. The file is removed
# when the R session ends.
hellinger_transform <- function(x) {
  x <- dplyr::mutate(x, quant = quant / sum(quant, na.rm = TRUE), .by = sample)
  path <- tempfile(fileext = ".parquet")
  duckplyr::compute_parquet(x, path)
  duckplyr::read_sql_duckdb(
    sprintf("SELECT sample, taxon, sqrt(quant) AS quant FROM read_parquet('%s')", gsub("'", "''", path)),
    prudence = "stingy"
  )
}

# duckplyr translates neither floor() nor a cross join, so "not a whole number"
# is round(x) != x and the table-wide values are window functions.
log_transform <- function(x, logbase) {
  x <- dplyr::mutate(
    x,
    minpos = min(quant, na.rm = TRUE),
    fractional = max(dplyr::if_else(round(quant) != quant, 1L, 0L), na.rm = TRUE)
  )
  x <- dplyr::mutate(x, quant = log(dplyr::if_else(fractional == 1L, quant / minpos, quant)) / log(logbase) + 1)
  dplyr::select(x, sample, taxon, quant)
}

check_quantities <- function(x, quant) {
  bad <- as.data.frame(dplyr::summarise(
    x,
    missing = sum(dplyr::if_else(is.na(quant), 1L, 0L)),
    negative = sum(dplyr::if_else(quant < 0, 1L, 0L), na.rm = TRUE)
  ))
  if (isTRUE(bad$missing > 0)) stop("`", quant, "` has missing values", call. = FALSE)
  if (isTRUE(bad$negative > 0)) stop("`", quant, "` has negative values", call. = FALSE)
  invisible()
}
