data(varespec, package = "vegan")

fixtures <- list(varespec = as.matrix(varespec), synthetic = synthetic_community(), lmo = lmo_subset()$matrix)
methods <- c("manhattan", "canberra")

for (fixture in names(fixtures)) {
  for (method in methods) {
    for (binary in c(FALSE, TRUE)) {
      test_that(paste0(method, ", binary = ", binary, " on ", fixture, " matches vegan::vegdist"), {
        m <- fixtures[[fixture]]
        got <- el_vegdist(long_table(m), method = method, binary = binary)

        expect_equal(as_full_matrix(got), as_full_matrix(vegan_reference(m, method, binary)), tolerance = 1e-10)
      })
    }
  }
}

for (method in methods) {
  for (binary in c(FALSE, TRUE)) {
    test_that(paste0(method, ", binary = ", binary, ": two all-zero samples listed in `samples` match vegan"), {
      m <- synthetic_community()
      empty <- c("empty1", "empty2")
      got <- el_vegdist(long_table(m), method = method, binary = binary, samples = c(rownames(m), empty))

      expect_equal(as_full_matrix(got), as_full_matrix(vegan_reference(m, method, binary, empty)), tolerance = 1e-10)
    })
  }
}

test_that("zero samples: manhattan to the other sample is its total and zero-zero is 0", {
  m <- synthetic_community()[1:3, ]
  got <- as.matrix(el_vegdist(long_table(m), method = "manhattan", samples = c(rownames(m), "e1", "e2")))

  expect_equal(unname(got["e1", rownames(m)]), unname(rowSums(m)))
  expect_equal(got["e1", "e2"], 0)
})

test_that("zero samples: canberra to another sample is 1 and zero-zero is NA, as in vegan", {
  m <- synthetic_community()[1:3, ]
  got <- as.matrix(el_vegdist(long_table(m), method = "canberra", samples = c(rownames(m), "e1", "e2")))

  expect_equal(unname(got["e1", rownames(m)]), c(1, 1, 1))
  expect_true(is.na(got["e1", "e2"]))
})

test_that("explicit zero rows are ignored for manhattan and canberra", {
  m <- synthetic_community()
  with_zeros <- rbind(long_table(m), data.frame(sample = "s04", taxon = "t59", quant = 0))

  for (method in methods) {
    for (binary in c(FALSE, TRUE)) {
      expect_equal(
        as_full_matrix(el_vegdist(with_zeros, method = method, binary = binary)),
        as_full_matrix(el_vegdist(long_table(m), method = method, binary = binary))
      )
    }
  }
})

test_that("negative quantities are an error for manhattan and canberra", {
  x <- long_table(synthetic_community())

  for (method in methods) expect_error(el_vegdist(transform(x, quant = -quant), method = method), "negative")
})

test_that("the result does not depend on batch_rows or threads", {
  x <- long_table(fixtures$lmo)

  for (method in methods) {
    reference <- el_vegdist(x, method = method, threads = 1)
    expect_equal(as.numeric(el_vegdist(x, method = method, batch_rows = 1)), as.numeric(reference), tolerance = 1e-12)
    expect_equal(as.numeric(el_vegdist(x, method = method, threads = 4)), as.numeric(reference), tolerance = 1e-12)
  }
})

test_that("integer counts above the 32-bit limit give vegan's result for manhattan and canberra", {
  set.seed(3)
  m <- matrix(0L, 10, 15, dimnames = list(sprintf("s%02d", 1:10), sprintf("t%02d", 1:15)))
  m[sample(length(m), 60)] <- sample(c(60000L, 3L, 2000000000L), 60, TRUE)
  x <- long_table_nonzero(m)
  x$quant <- as.integer(x$quant)

  for (method in methods) {
    expect_equal(as_full_matrix(el_vegdist(x, method = method)), as_full_matrix(vegan::vegdist(m, method)), tolerance = 1e-10)
  }
})

naive_canberra <- function(sample, taxon, quant, n) {
  acc <- numeric(n * (n - 1) / 2)
  count <- numeric(n * (n - 1) / 2)
  for (t in unique(taxon)) {
    rows <- which(taxon == t)
    for (a in rows) for (b in rows) {
      if (sample[a] < sample[b]) {
        k <- sample[a] * (2 * n - sample[a] - 1) / 2 + (sample[b] - sample[a] - 1) + 1
        acc[k] <- acc[k] + abs(quant[a] - quant[b]) / (quant[a] + quant[b])
        count[k] <- count[k] + 1
      }
    }
  }
  list(acc = acc, count = count)
}

test_that("accumulate_pairs in canberra mode adds |a - b| / (a + b) and counts the shared taxa", {
  d <- random_community(n_samples = 25L, n_taxa = 60L, seed = 17)
  expected <- naive_canberra(d$sample, d$taxon, d$quant, 25L)

  for (threads in c(1L, 3L)) {
    acc <- numeric(25 * 24 / 2)
    count <- numeric(25 * 24 / 2)
    accumulate_pairs(acc, d$sample, d$taxon, d$quant, 25L, canberra = TRUE, count = count, threads = threads)
    expect_equal(acc, expected$acc)
    expect_equal(count, expected$count)
  }
})
