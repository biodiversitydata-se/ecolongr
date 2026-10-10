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

synthetic_community <- function() {
  set.seed(42)
  m <- matrix(0, nrow = 20, ncol = 60, dimnames = list(sprintf("s%02d", 1:20), sprintf("t%02d", 1:60)))
  m[sample(length(m), 250)] <- rpois(250, 4) + 1
  m[2, ] <- m[1, ]
  m[3, ] <- 0
  m[3, "t05"] <- 7
  m[, "t60"] <- 0
  m[9, "t60"] <- 3
  m[rowSums(m) == 0, "t01"] <- 1
  m
}

vegan_reference <- function(m, method, binary, empty = character()) {
  full <- rbind(m, matrix(0, length(empty), ncol(m), dimnames = list(empty, colnames(m))))
  suppressWarnings(vegan::vegdist(full, method = method, binary = binary))
}
