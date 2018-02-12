#' Gather and bind the results of the benchmark jobs
#'
#' @param out_dir The folder in which to output intermediate and final results.
#' @param load_models Whether or not to load the models as well.
#'
#' @importFrom PRISM qsub_retrieve qacct qstat_j
#' @importFrom readr read_rds write_rds
#' @export
bs_bind_results <- function(out_dir, load_models = FALSE) {
  method_names <- list.dirs(out_dir, full.names = FALSE, recursive = FALSE) %>% discard(~ . == "")

  map_df(method_names, function(method_name) {
    method_folder <- paste0(out_dir, method_name)
    output_metrics_file <- paste0(method_folder, "/output.rds")
    output_models_file <- paste0(method_folder, "/output.rds")
    qsubhandle_file <- paste0(method_folder, "/qsubhandle.rds")

    if (file.exists(output_metrics_file)) {
      cat(method_name, ": Reading previous output\n", sep = "")
      output <- readr::read_rds(output_metrics_file)

      if (file.exists(output_models_file)) {
        models <- readr::read_rds(output_models_file)
        output$model <- models
      }

      output
    } else {
      cat(method_name, ": Output not found, skipping\n", sep = "")
      NULL
    }
  })
}
