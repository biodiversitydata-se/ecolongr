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

batch_table <- function() {
  data.frame(
    sample = c("a", "b", "c", "a", "b", "a", "c", "d", "e", "e"),
    taxon  = c("t3", "t3", "t3", "t1", "t1", "t2", "t2", "t2", "t2", "t4"),
    quant  = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  )
}

wide_from_long <- function(d, template, sample = "sample", taxon = "taxon", quant = "quant") {
  m <- matrix(0, nrow(template), ncol(template), dimnames = dimnames(template))
  m[cbind(match(d[[sample]], rownames(template)), match(d[[taxon]], colnames(template)))] <- d[[quant]]
  m
}

vegan_decostand <- function(m, method, ...) {
  out <- if (method == "rclr") {
    vegan::decostand(m, method = "rclr", impute = FALSE)
  } else {
    suppressWarnings(vegan::decostand(m, method = method, ...))
  }
  out <- as.matrix(out)
  attr(out, "decostand") <- NULL
  attr(out, "parameters") <- NULL
  out
}

long_table_nonzero <- function(m) {
  d <- data.frame(
    sample = rep(rownames(m), times = ncol(m)),
    taxon = rep(colnames(m), each = nrow(m)),
    quant = as.vector(as.matrix(m))
  )
  d[d$quant != 0, ]
}

align_signs <- function(a, b) {
  a <- as.matrix(a)
  b <- as.matrix(b)
  a * rep(sign(colSums(a * b)), each = nrow(a))
}

big_sparse_community <- function(n = 400, p = 1500, seed = 5) {
  set.seed(seed)
  m <- matrix(0, n, p, dimnames = list(sprintf("s%03d", seq_len(n)), sprintf("t%04d", seq_len(p))))
  k <- 0.04 * n * p
  idx <- sample(length(m), k)
  m[idx] <- rpois(k, 4) + 1
  m[, colSums(m) > 0, drop = FALSE]
}

compare_rda <- function(got, m, ..., k = 3, tolerance = 1e-8) {
  ref <- vegan::rda(m, ...)
  keep <- intersect(colnames(m), rownames(vegan::scores(got, choices = 1, display = "species")))

  expect_equal(unname(vegan::eigenvals(got)[seq_len(k)]), unname(vegan::eigenvals(ref)[seq_len(k)]), tolerance = tolerance)
  expect_equal(got$tot.chi, ref$tot.chi, tolerance = tolerance)
  for (scaling in 1:3) {
    s_got <- vegan::scores(got, choices = seq_len(k), scaling = scaling)
    s_ref <- vegan::scores(ref, choices = seq_len(k), scaling = scaling)
    sites_ref <- s_ref$sites[rownames(s_got$sites), , drop = FALSE]
    species_ref <- s_ref$species[keep, , drop = FALSE]
    expect_equal(unname(s_got$sites), unname(align_signs(sites_ref, s_got$sites)), tolerance = tolerance)
    expect_equal(unname(s_got$species[keep, , drop = FALSE]), unname(align_signs(species_ref, s_got$species[keep, , drop = FALSE])), tolerance = tolerance)
  }
}

lmo_subset <- function() {
  d <- utils::read.delim(gzfile(testthat::test_path("fixtures", "lmo_subset.tsv.gz")), stringsAsFactors = FALSE)
  long <- data.frame(sample = d$sample, taxon = d$asv_id, quant = d$count)
  m <- matrix(0, length(unique(long$sample)), length(unique(long$taxon)),
              dimnames = list(sort(unique(long$sample)), sort(unique(long$taxon))))
  m[cbind(match(long$sample, rownames(m)), match(long$taxon, colnames(m)))] <- long$quant
  list(long = long, matrix = m)
}

present_taxa <- function(m) m[, colSums(m) > 0, drop = FALSE]

clr_reference <- function(m, pseudocount) {
  pc <- if (length(pseudocount) == 1) rep(pseudocount, nrow(m)) else pseudocount[rownames(m)]
  y <- log(m + pc)
  y - rowMeans(y)
}

shuffled_within_taxon <- function(d, seed = 3) {
  set.seed(seed)
  d[order(d$taxon, runif(nrow(d))), ]
}

random_community <- function(n_samples = 60L, n_taxa = 300L, seed = 11) {
  set.seed(seed)
  d <- unique(data.frame(
    sample = sample.int(n_samples, 1500, TRUE) - 1L,
    taxon = sort(sample.int(n_taxa, 1500, TRUE))
  ))
  d$quant <- rpois(nrow(d), 5) + 1
  shuffled_within_taxon(d)
}

rclr_triplets <- function(m) {
  r <- suppressWarnings(vegan::decostand(m, "rclr", impute = FALSE))
  r <- as.matrix(r)
  idx <- which(!is.na(r), arr.ind = TRUE)
  list(i = idx[, 1] - 1L, j = idx[, 2] - 1L, v = r[idx], n = nrow(m), p = ncol(m))
}

double_centre <- function(m) {
  m <- sweep(m, 2, colMeans(m))
  sweep(m, 1, rowMeans(m))
}

vegan_rclr_imputed <- function(m) {
  out <- as.matrix(suppressWarnings(vegan::decostand(m, "rclr", impute = TRUE)))
  matrix(as.numeric(out), nrow(out), ncol(out), dimnames = dimnames(out))
}
