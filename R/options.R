resolve_threads <- function(threads) {
  if (is.null(threads)) threads <- getOption("ecolongr.threads", max(1L, parallel::detectCores() %/% 2L))
  if (!is.numeric(threads) || length(threads) != 1 || is.na(threads) || threads < 1 || threads != round(threads)) {
    stop("`threads` must be a positive whole number", call. = FALSE)
  }
  as.integer(threads)
}

resolve_memory_limit <- function(memory_limit) {
  if (is.null(memory_limit)) memory_limit <- getOption("ecolongr.memory_limit", "250MB")
  if (!is.character(memory_limit) || length(memory_limit) != 1 ||
      !grepl("^[0-9]+(\\.[0-9]+)? ?[KMGT]i?B$", memory_limit)) {
    stop("`memory_limit` must be a single size such as \"250MB\"", call. = FALSE)
  }
  memory_limit
}
