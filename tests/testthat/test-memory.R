# Peak memory of el_vegdist() on a skewed table, measured in a fresh R process.
# Pairs per taxon grow with the square of the samples that share it, so a
# few common taxa make the number of sample pairs large.
peak_extra_mb <- function(n_samples, n_taxa, n_draws) {
  pkg <- system.file(package = "ecolongr")
  script <- tempfile(fileext = ".R")
  writeLines(sprintf('
    suppressMessages(if (dir.exists(file.path("%s", "R"))) pkgload::load_all("%s", quiet = TRUE) else library(ecolongr))
    set.seed(1)
    p <- 1 / (1:%d)^1; p <- p / sum(p)
    d <- unique(data.frame(sample = sprintf("s%%05d", sample.int(%d, %d, TRUE)),
                           taxon = sprintf("t%%06d", sample.int(%d, %d, TRUE, prob = p))))
    d$quant <- rpois(nrow(d), 5) + 1
    path <- tempfile(fileext = ".parquet")
    con <- DBI::dbConnect(duckdb::duckdb()); DBI::dbWriteTable(con, "d", d)
    DBI::dbExecute(con, sprintf("COPY d TO \'%%s\' (FORMAT PARQUET)", path)); DBI::dbDisconnect(con, shutdown = TRUE)
    rm(d, p, con); invisible(gc())
    kb <- function(field) as.numeric(sub("[^0-9]*([0-9]+).*", "\\\\1", grep(paste0("^", field), readLines("/proc/self/status"), value = TRUE)))
    writeLines("5", "/proc/self/clear_refs")
    before <- kb("VmRSS")
    x <- duckplyr::read_parquet_duckdb(path)
    d <- el_vegdist(x)
    cat("EXTRA_MB", (kb("VmHWM") - before) / 1024, "\\n")
  ', pkg, pkg, n_taxa, n_samples, n_draws, n_taxa, n_draws), script)
  out <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), script, stdout = TRUE, stderr = FALSE))
  as.numeric(sub("EXTRA_MB ", "", grep("^EXTRA_MB", out, value = TRUE)))
}

test_that("el_vegdist on a skewed table with 73M sample pairs adds under 400 MB", {
  skip_on_cran()
  skip_if_not(file.exists("/proc/self/clear_refs"), "needs Linux /proc")
  skip_if_not_installed("pkgload")

  extra <- peak_extra_mb(n_samples = 1000, n_taxa = 50000, n_draws = 600000)

  expect_length(extra, 1)
  expect_lt(extra, 400)
})
