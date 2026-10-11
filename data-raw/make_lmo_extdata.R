# Builds inst/extdata/lmo_free_living*.tsv* from the published LMO 16S dataset.
# Download https://www.gbif.se/ipt/archive.do?r=lmo_2011-2023_16s, unzip it, and run:
#   Rscript data-raw/make_lmo_extdata.R path/to/unzipped/archive
library(duckplyr)

archive <- commandArgs(TRUE)[1]
tsv <- list(delim = "\t", quote = "", header = TRUE)
occ <- read_csv_duckdb(file.path(archive, "occurrence.txt"), options = tsv, prudence = "stingy")
emof <- read_csv_duckdb(file.path(archive, "extendedmeasurementorfact.txt"),
                        options = c(tsv, list(all_varchar = TRUE)), prudence = "stingy")

fraction <- emof |>
  filter(measurementType == "filter_fraction", measurementValue == "3-0.2") |>
  select(occurrenceID)

counts <- occ |>
  semi_join(fraction, by = "occurrenceID") |>
  transmute(sample = eventID, taxon = taxonID, quant = organismQuantity)

n_samples <- counts |> summarise(n = n_distinct(sample)) |> pull(n)
common <- counts |>
  summarise(prevalence = n_distinct(sample), .by = taxon) |>
  filter(prevalence >= 0.05 * n_samples) |>
  select(taxon)

kept <- counts |> semi_join(common, by = "taxon") |> as.data.frame()
kept$taxon <- sub("^ASV:", "", kept$taxon)
kept <- kept[order(kept$sample, kept$taxon), ]

samples <- occ |>
  semi_join(fraction, by = "occurrenceID") |>
  distinct(sample = eventID, date = eventDate, depth = minimumDepthInMeters) |>
  as.data.frame()

measured <- emof |>
  semi_join(fraction, by = "occurrenceID") |>
  filter(measurementType %in% c("temperature", "salinity")) |>
  inner_join(select(occ, occurrenceID, sample = eventID), by = "occurrenceID") |>
  distinct(sample, measurementType, measurementValue) |>
  as.data.frame()
for (type in c("temperature", "salinity")) {
  one <- measured[measured$measurementType == type, ]
  samples[[type]] <- as.numeric(one$measurementValue[match(samples$sample, one$sample)])
}

samples <- samples[samples$sample %in% kept$sample, ]
samples <- samples[order(samples$date, samples$sample), ]

con <- gzfile("inst/extdata/lmo_free_living.tsv.gz", "w", compression = 9)
utils::write.table(kept, con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)
utils::write.table(samples, "inst/extdata/lmo_free_living_samples.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
