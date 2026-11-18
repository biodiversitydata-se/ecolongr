library(dplyr)

# A long table of read counts for 16S ASVs in a time series from the Baltic Sea.
counts <- duckplyr::read_csv_duckdb(
  system.file("extdata", "lmo_free_living.tsv.gz", package = "ecolongr"),
  options = list(delim = "\t", header = TRUE), prudence = "stingy"
)
metadata <- read.delim(system.file("extdata", "lmo_free_living_samples.tsv", package = "ecolongr"))

# Samples with at least 5,000 reads.
deep <- counts |>
  summarise(reads = sum(quant, na.rm = TRUE), .by = sample) |>
  filter(reads >= 5000) |>
  select(sample)
counts <- semi_join(counts, deep, by = "sample")
metadata <- filter(metadata, sample %in% as.data.frame(deep)$sample)
