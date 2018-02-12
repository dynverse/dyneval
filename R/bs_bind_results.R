#' Gather and bind the results of the benchmark jobs
#'
#' @param out_dir The folder in which to output intermediate and final results.
#' @param load_models Whether or not to load the models as well.
#'
#' @importFrom readr read_rds
#' @export
bs_bind_results <- function(out_dir, load_models = FALSE) {
  method_names <- list.dirs(out_dir, full.names = FALSE, recursive = FALSE) %>% discard(~ . == "")

  # process each method separately
  map_df(method_names, function(method_name) {
    method_folder <- paste0(out_dir, method_name)
    output_metrics_file <- paste0(method_folder, "/output_metrics.rds")
    output_models_file <- paste0(method_folder, "/output_models.rds")

    if (file.exists(output_metrics_file)) {
      cat(method_name, ": Reading previous output\n", sep = "")
      output <- readr::read_rds(output_metrics_file)

      # read models, if requested
      if (file.exists(output_models_file)) {
        models <- readr::read_rds(output_models_file)
        output$model <- models
      }

      output
    } else {
      cat(method_name, ": Output not found, skipping\n", sep = "")
      NULL
    }
  }) %>%
    as_tibble()
}
