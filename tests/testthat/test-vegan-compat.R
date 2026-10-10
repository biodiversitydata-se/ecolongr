data(varespec, varechem, package = "vegan")

m <- as.matrix(varespec)
group <- cut(varechem$pH, 3, labels = c("low", "mid", "high"))
env <- varechem[, c("N", "Humdepth")]
metadata <- data.frame(group = group, row.names = rownames(varechem))

with_pdf <- function(code) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  force(code)
}

dist_pairs <- list(
  bray = list(
    got = function() el_vegdist(long_table(m), "bray", samples = rownames(m)),
    ref = function() vegan::vegdist(m, "bray")
  ),
  aitchison = list(
    got = function() el_vegdist(long_table(m), "aitchison", samples = rownames(m)),
    ref = function() vegan::vegdist(m, "aitchison", pseudocount = 1)
  )
)

for (method in names(dist_pairs)) {
  test_that(paste0("a ", method, " dist behaves as vegan's in metaMDS, adonis2, betadisper, hclust and cmdscale"), {
    got <- dist_pairs[[method]]$got()
    ref <- dist_pairs[[method]]$ref()

    set.seed(1)
    nmds_got <- vegan::metaMDS(got, k = 2, trace = 0, try = 3, trymax = 3)
    set.seed(1)
    nmds_ref <- vegan::metaMDS(ref, k = 2, trace = 0, try = 3, trymax = 3)
    expect_s3_class(nmds_got, "metaMDS")
    expect_equal(nmds_got$stress, nmds_ref$stress, tolerance = 1e-6)
    expect_lt(vegan::procrustes(nmds_ref, nmds_got)$ss, 1e-6)

    set.seed(2)
    ad_got <- vegan::adonis2(got ~ group, data = metadata, permutations = 199)
    set.seed(2)
    ad_ref <- vegan::adonis2(ref ~ group, data = metadata, permutations = 199)
    expect_equal(ad_got$R2, ad_ref$R2, tolerance = 1e-8)
    expect_equal(ad_got$F, ad_ref$F, tolerance = 1e-8)
    expect_equal(ad_got$`Pr(>F)`, ad_ref$`Pr(>F)`)
    expect_lt(ad_got$`Pr(>F)`[1], 0.2)

    bd_got <- vegan::betadisper(got, group)
    bd_ref <- vegan::betadisper(ref, group)
    expect_equal(bd_got$distances, bd_ref$distances, tolerance = 1e-8)

    expect_equal(stats::hclust(got, "average")$height, stats::hclust(ref, "average")$height, tolerance = 1e-8)
    expect_equal(stats::cmdscale(got, k = 2, eig = TRUE)$eig, stats::cmdscale(ref, k = 2, eig = TRUE)$eig, tolerance = 1e-8)
  })
}

test_that("a dist from el_vegdist has the labels, size and attributes vegan's has, in the order of `samples`", {
  got <- dist_pairs$bray$got()
  ref <- dist_pairs$bray$ref()

  expect_equal(attr(got, "Labels"), rownames(m))
  for (a in c("Size", "Labels", "Diag", "Upper")) expect_equal(attr(got, a), attr(ref, a))
  expect_equal(as.matrix(got), as.matrix(ref), tolerance = 1e-10)
})

test_that("without `samples` the samples are in sorted order, and `samples` puts an ordination in the order of the metadata", {
  default <- el_rda(long_table(m), k = 2, solver = "exact")
  ordered <- el_rda(long_table(m), k = 2, solver = "exact", samples = rownames(m))

  expect_equal(rownames(vegan::scores(default, display = "sites")), sort(rownames(m)))
  expect_equal(rownames(vegan::scores(ordered, display = "sites")), rownames(m))
})

ordinations <- list(
  none = list(got = function() el_rda(long_table(m), k = 5, solver = "exact", samples = rownames(m)), ref = function() vegan::rda(m)),
  clr = list(
    got = function() el_rda(long_table(m), k = 5, transform = "clr", pseudocount = 1, solver = "exact", samples = rownames(m)),
    ref = function() vegan::rda(vegan::decostand(m, "clr", pseudocount = 1))
  ),
  rclr = list(
    got = function() el_rda(long_table(m), k = 3, transform = "rclr", samples = rownames(m)),
    ref = function() vegan::rda(vegan::decostand(m, "rclr"))
  )
)

for (transform in names(ordinations)) {
  test_that(paste0("an el_rda result (", transform, ") works with summary, plot, biplot, ordiplot, ordihull and ordiellipse"), {
    got <- ordinations[[transform]]$got()
    ref <- ordinations[[transform]]$ref()
    k <- length(vegan::eigenvals(got))

    importance <- summary(got)$cont$importance
    expect_equal(unname(importance[2, seq_len(k)]), unname(summary(ref)$cont$importance[2, seq_len(k)]), tolerance = 1e-8)

    with_pdf({
      expect_no_error(plot(got))
      expect_no_error(plot(got, display = "sites", choices = c(1, 3)))
      expect_no_error(suppressWarnings(stats::biplot(got)))
      ord <- vegan::ordiplot(got, display = "sites")
      expect_no_error(vegan::ordihull(ord, group))
      expect_no_error(vegan::ordiellipse(ord, group))
      expect_no_error(vegan::ordispider(ord, group))
    })
  })
}

test_that("envfit on an el_rda result gives vegan's r2 and p-values on the first two axes", {
  got <- el_rda(long_table(m), k = 5, solver = "exact", samples = rownames(m))
  ref <- vegan::rda(m)

  set.seed(3)
  fit_got <- vegan::envfit(got, env, choices = 1:2, permutations = 199)
  set.seed(3)
  fit_ref <- vegan::envfit(ref, env, choices = 1:2, permutations = 199)

  expect_equal(fit_got$vectors$r, fit_ref$vectors$r, tolerance = 1e-8)
  expect_equal(fit_got$vectors$pvals, fit_ref$vectors$pvals)
  expect_equal(abs(unname(fit_got$vectors$arrows)), abs(unname(fit_ref$vectors$arrows)), tolerance = 1e-8)
})

test_that("procrustes between an el_rda result and vegan's rda is zero, so the sites are the same configuration", {
  for (transform in c("none", "clr")) {
    got <- ordinations[[transform]]$got()
    ref <- ordinations[[transform]]$ref()

    expect_lt(vegan::procrustes(ref, got)$ss, 1e-8)
  }
})

test_that("functions that need the full matrix fail with an error, not a wrong answer", {
  got <- el_rda(long_table(m), k = 3, solver = "exact")

  expect_error(stats::predict(got))
  expect_error(stats::fitted(got))
})

test_that("anova of an unconstrained el_rda result gives vegan's table: no constrained component", {
  got <- el_rda(long_table(m), k = 3, solver = "exact", samples = rownames(m))
  ref <- vegan::rda(m)

  expect_equal(as.data.frame(stats::anova(got)), as.data.frame(stats::anova(ref)), ignore_attr = TRUE, tolerance = 1e-8)
})
