filename <- paste0(tempdir(), "/dyneval_test_datasets.rds")
.datasets_location = "../../../dyngen/results/4/"
set.seed(1)
datasets_sel <- dplyr::sample_n(load_datasets_info(), 16)
cat("Loading ", nrow(datasets_sel), " datasets\n", sep="")
datasets <- load_datasets(mc_cores = 8, datasets_sel)
saveRDS(datasets, filename)
