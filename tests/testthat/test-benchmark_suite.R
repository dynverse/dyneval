# context("Benchmark suite")
#
# test_that(paste0("Check benchmark suite"), {
#   tmpdir <- tempfile("benchmarksuitetest")
#   dir.create(tmpdir)
#
#   tasks <- dyntoy::generate_toy_datasets(ti_types = "bifurcating", num_replicates = 2, num_cells = 30, num_genes = 31)
#
#   out <- benchmark_suite_submit(
#     tasks = tasks,
#     task_group = tasks$type,
#     task_fold = tasks$replicate,
#     out_dir = tmpdir,
#     timeout = 120,
#     methods = get_descriptions(as_tibble = TRUE) %>% filter(short_name == "CTVEM"),
#     metrics = c("auc_R_nx"),
#     num_cores = 4,
#     num_iterations = 2,
#     num_init_params = 10,
#     num_repeats = 1,
#     do_it_local = TRUE
#   )
#
#   unlink(tmpdir, recursive = TRUE)
#   # # # expect to fail
#   # out <- execute_evaluation(tasks = dyntoy::toy_tasks[5,], method = description_celltree_maptpx(), parameters = list(), timeout = 4, metrics = "auc_R_nx")
#   # attr(out, "extras")$.summary
# })
#
# test_that(paste0("Check benchmark suite"), {
#   tmpdir <- tempfile("benchmarksuitetest")
#   dir.create(tmpdir)
#
#   tasks <- dyntoy::generate_toy_datasets(ti_types = "bifurcating", num_replicates = 2, num_cells = 30, num_genes = 31)
#
#   out <- benchmark_suite_submit(
#     tasks = tasks,
#     task_group = tasks$type,
#     task_fold = tasks$replicate,
#     out_dir = tmpdir,
#     timeout = 0,
#     methods = get_descriptions(as_tibble = TRUE) %>% filter(short_name == "CTVEM"),
#     metrics = c("auc_R_nx"),
#     num_cores = 4,
#     num_iterations = 2,
#     num_init_params = 10,
#     num_repeats = 1,
#     do_it_local = TRUE
#   )
#
#   unlink(tmpdir, recursive = TRUE)
#   # # # expect to fail
#   # out <- execute_evaluation(tasks = dyntoy::toy_tasks[5,], method = description_celltree_maptpx(), parameters = list(), timeout = 4, metrics = "auc_R_nx")
#   # attr(out, "extras")$.summary
# })
