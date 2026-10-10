# usage: Rscript bench/run.R table.parquet threads [memory_limit]
# Use an installed package (R CMD INSTALL), not devtools::load_all(), which compiles without optimisation.
library(ecolongr)
a <- commandArgs(TRUE)
x <- duckplyr::read_parquet_duckdb(a[1])
limit <- if (length(a) >= 3) a[3] else NULL
t <- system.time(d <- el_vegdist(x, threads = as.integer(a[2]), memory_limit = limit))
peak <- as.numeric(sub("[^0-9]*([0-9]+).*", "\\1", grep("^VmHWM", readLines("/proc/self/status"), value = TRUE))) / 1024
cat(sprintf("n=%d threads=%s elapsed %.1fs peak RSS %.0f MB\n", attr(d, "Size"), a[2], t[["elapsed"]], peak))
