# DuckDB splits a table of more than one row group (122,880 rows) over its threads,
# so these tests need a table that large; the small fixtures never exercise that.

test_that("the sorted table is in taxon order when DuckDB uses several threads", {
  x <- many_taxa_table()
  expect_gt(nrow(x), 122880)
  samples <- sort(unique(x$sample))

  for (run in 1:5) {
    stream <- open_taxon_stream(x, threads = 4)
    taxa <- integer()
    stream_batches(stream, samples, 1e9, function(b) taxa <<- c(taxa, b$taxon))
    close_taxon_stream(stream)

    expect_false(is.unsorted(taxa))
    expect_identical(length(taxa), nrow(x))
  }
})

test_that("el_vegdist matches vegan with several threads, run after run", {
  x <- many_taxa_table()
  m <- wide_from_long(x, matrix(0, 95, 4500, dimnames = list(sprintf("s%03d", 1:95), sprintf("t%04d", 1:4500))))
  reference <- list(
    bray = vegan::vegdist(m, "bray"),
    jaccard = vegan::vegdist(m, "jaccard", binary = TRUE),
    euclidean = vegan::vegdist(m, "euclidean"),
    manhattan = vegan::vegdist(m, "manhattan"),
    canberra = suppressWarnings(vegan::vegdist(m, "canberra")),
    aitchison = stats::dist(clr_reference(m, 1))
  )

  for (run in 1:3) {
    for (method in names(reference)) {
      got <- el_vegdist(x, method = method, binary = method == "jaccard", samples = rownames(m), threads = 4)
      expect_equal(as.numeric(got), as.numeric(reference[[method]]), tolerance = 1e-10,
                   info = paste(method, "run", run))
    }
  }
})

test_that("el_rda matches vegan on a table large enough for parallel DuckDB", {
  x <- many_taxa_table()
  m <- wide_from_long(x, matrix(0, 95, 4500, dimnames = list(sprintf("s%03d", 1:95), sprintf("t%04d", 1:4500))))

  compare_rda(el_rda(x, k = 3, samples = rownames(m)), m, k = 3)
})
