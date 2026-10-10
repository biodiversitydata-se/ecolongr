# Builds lmo_subset.tsv.gz and lmo_subset_samples.tsv from the LMO 16S time series.
# usage: Rscript make_lmo_subset.R lmo16s.tsv.gz lmo_samples.tsv
suppressMessages(library(dplyr))
a <- commandArgs(TRUE)
samples <- read.delim(a[2], stringsAsFactors = FALSE, na.strings = c("NA", "")) |>
  filter(project == "LMO", seqtype == "16SDNA", filter_fraction == 0.2, !replaced)
counts <- read.delim(gzfile(a[1]), stringsAsFactors = FALSE) |> filter(sample %in% samples$sample)
samples <- samples |> filter(sample %in% counts$sample) |> arrange(date)
keep <- samples[round(seq(1, nrow(samples), length.out = 90)), ]
sub <- counts |> filter(sample %in% keep$sample)
asvs <- sub |> count(asv_id) |> filter(n >= 8) |> pull(asv_id)
sub <- sub |> filter(asv_id %in% asvs) |> arrange(asv_id, sample)
write.table(sub, gzfile("lmo_subset.tsv.gz"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(keep |> select(sample, date, depth_m), "lmo_subset_samples.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat(nrow(sub), "rows,", n_distinct(sub$sample), "samples,", n_distinct(sub$asv_id), "ASVs\n")
