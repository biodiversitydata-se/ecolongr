long_table <- function(m) {
  d <- data.frame(
    sample = rep(rownames(m), times = ncol(m)),
    taxon = rep(colnames(m), each = nrow(m)),
    quant = as.vector(as.matrix(m))
  )
  d[d$quant > 0, ]
}

as_full_matrix <- function(d) {
  m <- as.matrix(d)
  m[order(rownames(m)), order(colnames(m))]
}
