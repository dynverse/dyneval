# where are the dyngen datasets stored?
.datasets_location = "../../../dyngen/results/4/"
# .datasets_location = "../dyngen/results/4/"

# set a seed for replication reasons
set.seed(1)

# sample some datasets
task_sel <- dplyr::sample_n(load_datasets_info(), 8)

# load them
cat("Loading ", nrow(task_sel), " datasets\n", sep="")
tasks <- load_datasets(mc_cores = 8, task_sel)

# save file
saveRDS(tasks, paste0(tempdir(), "/dyneval_test_datasets.rds"))
