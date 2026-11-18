utils::globalVariables(c("h"))

#' Alpha diversity of each sample
#'
#' Counterparts of [vegan::diversity()] and [vegan::specnumber()] for a long
#' table without zeros. Each is one grouped aggregate in DuckDB, and only the
#' result, one number per sample, is read into R.
#'
#' * `el_diversity()`: the Shannon index `-sum(p * log(p))` (divided by
#'   `log(base)`), the Simpson index `1 - sum(p * p)` or the inverse Simpson
#'   index `1 / sum(p * p)`, where `p` is the quantity divided by the sample
#'   total. As in `vegan`, quantities need not be whole numbers, and a common
#'   scale does not change the result, so proportions are accepted.
#' * `el_specnumber()`: the number of taxa with a quantity above zero.
#'
#' Rows with `quant == 0` are dropped. A sample with no non-zero row cannot be
#' told from a sample that is not in the table, so it is absent from the result.
#' `vegan` returns for an all-zero row Shannon 0, Simpson 1, inverse Simpson
#' `Inf` and richness 0; add these entries yourself if you need them.
#'
#' Not implemented: `groups`, `MARGIN` and `equalize.groups` of `diversity()`,
#' and `renyi()`, `fisher.alpha()`, `estimateR()`, `specpool()` and evenness.
#'
#' @inheritParams el_decostand
#' @param index One of `"shannon"`, `"simpson"` and `"invsimpson"`.
#' @param base Base of the logarithm for `"shannon"`; ignored by the others.
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
  x <- dplyr::mutate(x, quant = quant / sum(quant, na.rm = TRUE), .by = sample)
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
  stats::setNames(values, sample)[order(sample)]
}
