# usage: Rscript bench/gen.R n_samples n_taxa n_draws zipf_exponent out.parquet
# Synthetic skewed table: a few taxa occur in nearly every sample, which is the expensive case.
a <- commandArgs(TRUE); ns <- as.integer(a[1]); nt <- as.integer(a[2]); nnz <- as.numeric(a[3]); z <- as.numeric(a[4])
set.seed(1)
p <- 1/(1:nt)^z; p <- p/sum(p)
# draw nnz samples/taxa pairs; duplicates removed
i <- sample.int(ns, nnz, TRUE); j <- sample.int(nt, nnz, TRUE, prob = p)
d <- unique(data.frame(sample = sprintf("s%05d", i), taxon = sprintf("t%06d", j)))
d$quant <- rpois(nrow(d), 5) + 1
cat("rows", nrow(d), " max samples/taxon", max(table(d$taxon)), " sum k^2 (M)", sum(as.numeric(table(d$taxon))^2)/1e6, "\n")
con <- DBI::dbConnect(duckdb::duckdb()); DBI::dbWriteTable(con, "d", d)
DBI::dbExecute(con, sprintf("COPY d TO '%s' (FORMAT PARQUET)", a[5])); DBI::dbDisconnect(con, shutdown = TRUE)
