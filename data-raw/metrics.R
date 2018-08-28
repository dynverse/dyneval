library(tibble)

metrics <- readr::read_tsv("data-raw/metrics.tsv", col_types = cols(perfect = "d", worst = "d", .default = "c"))

devtools::use_data(metrics, overwrite = TRUE)
