#' Implicit centred clr matrix
#'
#' The centred clr matrix is `B - u 1' - 1 c'` with `B` sparse, so `irlba` only
#' needs its products with vectors. `scale` holds the divisor of each column.
#'
#' @name clr_operator-class
#' @aliases clr_operator-class %*%,clr_operator,numeric-method %*%,numeric,clr_operator-method
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom methods setClass setMethod new as
#' @keywords internal
setClass("clr_operator", contains = "dgCMatrix", slots = c(u = "numeric", cvec = "numeric", scale = "numeric"))

setMethod("%*%", signature(x = "clr_operator", y = "numeric"), function(x, y) {
  y <- y / x@scale
  b <- as(x, "dgCMatrix")
  matrix(as.numeric(b %*% y) - x@u * sum(y) - sum(x@cvec * y), ncol = 1)
})

setMethod("%*%", signature(x = "numeric", y = "clr_operator"), function(x, y) {
  b <- as(y, "dgCMatrix")
  matrix((as.numeric(x %*% b) - sum(y@u * x) - sum(x) * y@cvec) / y@scale, nrow = 1)
})
