utils::globalVariables(c("h", "total"))

#' Alpha diversity of each sample
#'
#' Counterparts of [vegan::diversity()] and [vegan::specnumber()] for a long
#' table without zeros. The query runs in DuckDB and only the result, one number
#' per sample, is read into R; no taxa by samples matrix is built.
#' `el_diversity()` computes the sample totals, joins each total back onto its
#' rows and takes one grouped sum. `el_specnumber()` is one grouped count.
#'
#' * `el_diversity()`: the Shannon index `-sum(p * log(p))` (divided by
#'   `log(base)`), the Simpson index `1 - sum(p * p)` or the inverse Simpson
#'   index `1 / sum(p * p)`, where `p` is the quantity divided by the sample
#'   total. As in `vegan`, quantities need not be whole numbers, and a common
#'   scale does not change the result, so proportions are accepted.
#' * `el_specnumber()`: the number of rows with a quantity above zero in each
#'   sample.
#'
#' The table must have one row per sample and taxon. Rows are not merged, so a
#' sample and taxon that occur twice count as two taxa.
#'
#' Rows with `quant == 0` are dropped. A sample with no non-zero row cannot be
#' told from a sample that is not in the table, so it is absent from the result.
#' `vegan` returns for an all-zero row Shannon 0, Simpson 1, inverse Simpson
#' `Inf` and richness 0; add these entries yourself if you need them.
#'
#' The names of the result are the samples as character, sorted as text: factor
#' samples are not in level order, numeric samples sort as `"10" < "9"`, and an
#' `NA` sample is named `NA` and comes last.
#'
#' The sums are not compensated, so a sample with millions of taxa can differ
#' from `vegan` in the last digits (about 2e-11 relative with 1e7 equally common
#' taxa), and which digits can vary between runs with several threads.
#'
#' Not implemented: `groups`, `MARGIN` and `equalize.groups` of `diversity()`,
#' and `renyi()`, `fisher.alpha()`, `estimateR()`, `specpool()` and evenness.
#'
#' @param x A data frame or lazy `duckplyr` table with one row per non-zero
#'   sample and taxon. Rows with `quant == 0` are dropped.
#' @param index One of `"shannon"`, `"simpson"` and `"invsimpson"`.
#' @param base Base of the logarithm for `"shannon"`; ignored by the others.
#' @param check If `TRUE`, run one query to stop on negative, missing or
#'   non-finite (NaN, Inf) quantities. A sample total above about 1e308 is not
#'   checked. Use `FALSE` to skip the query.
#' @param sample,taxon,quant Names of the columns in `x`.
#' @return A numeric vector with one value per sample, named by sample and
#'   sorted by name.
#' @export
el_diversity <- function(x, index = c("shannon", "simpson", "invsimpson"), base = exp(1),
                         check = TRUE, sample = "sample", taxon = "taxon", quant = "quant") {
  index <- match.arg(index)
  if (!is.numeric(base) || length(base) != 1 || !is.finite(base) || base <= 0 || base == 1) {
    stop("`base` must be a single positive number other than 1", call. = FALSE)
  }
  x <- alpha_input(x, check, sample, taxon, quant)
  totals <- dplyr::summarise(x, total = sum(quant, na.rm = TRUE), .by = sample)
  x <- dplyr::inner_join(x, totals, by = "sample")
  x <- dplyr::mutate(x, quant = quant / total)
  x <- if (index == "shannon") {
    dplyr::summarise(x, h = sum(quant * log(quant), na.rm = TRUE), .by = sample)
  } else {
    dplyr::summarise(x, h = sum(quant * quant, na.rm = TRUE), .by = sample)
  }
  d <- as.data.frame(x)
  h <- switch(index, shannon = -d$h / log(base), simpson = 1 - d$h, invsimpson = 1 / d$h)
  by_sample(h, d$sample)
}

#' @rdname el_diversity
#' @export
el_specnumber <- function(x, check = TRUE, sample = "sample", taxon = "taxon", quant = "quant") {
  x <- alpha_input(x, check, sample, taxon, quant)
  d <- as.data.frame(dplyr::summarise(x, k = dplyr::n(), .by = sample))
  by_sample(d$k, d$sample)
}

alpha_input <- function(x, check, sample, taxon, quant) {
  x <- duckplyr::as_duckdb_tibble(x, prudence = "stingy")
  x <- dplyr::select(x, sample = dplyr::all_of(sample), taxon = dplyr::all_of(taxon),
                     quant = dplyr::all_of(quant))
  if (check) check_quantities(x, quant)
  dplyr::filter(x, quant > 0)
}

by_sample <- function(values, sample) {
  sample <- as.character(sample)
  stats::setNames(values, sample)[order(sample)]
}
