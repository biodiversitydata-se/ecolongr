check_pseudocount_form <- function(pseudocount) {
  inverse_depth <- is.character(pseudocount) && length(pseudocount) == 1 && identical(pseudocount, "inverse_depth")
  scalar <- is.numeric(pseudocount) && length(pseudocount) == 1 && is.null(names(pseudocount))
  named <- is.numeric(pseudocount) && length(pseudocount) >= 1 && !is.null(names(pseudocount)) && !anyNA(names(pseudocount))
  ok <- inverse_depth || ((scalar || named) && all(is.finite(pseudocount)) && all(pseudocount > 0))
  if (!ok) {
    stop("`pseudocount` must be a positive number, \"inverse_depth\", or a positive numeric vector named by sample",
         call. = FALSE)
  }
  invisible()
}

# One pseudocount per sample, in the units of the counts. "inverse_depth" is 1 / the
# sample's total count; an all-zero sample gets 1, which does not matter for its clr.
resolve_pseudocount <- function(pseudocount, totals) {
  check_pseudocount_form(pseudocount)
  samples <- as.character(totals$sample)
  if (identical(pseudocount, "inverse_depth")) {
    pc <- ifelse(totals$total > 0, 1 / totals$total, 1)
  } else if (is.null(names(pseudocount))) {
    pc <- rep(pseudocount, length(samples))
  } else {
    missing <- setdiff(samples, names(pseudocount))
    if (length(missing) > 0) {
      stop("`pseudocount` has no pseudocount for sample(s): ", paste(utils::head(missing, 5), collapse = ", "), call. = FALSE)
    }
    pc <- unname(pseudocount[samples])
  }
  stats::setNames(pc, samples)
}
