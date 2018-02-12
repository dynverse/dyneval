#' Fetch the results of the benchmark jobs from the cluster.
#'
#' @param out_dir The folder in which to output intermediate and final results.
#'
#' @importFrom PRISM qsub_retrieve qacct qstat_j
#' @importFrom readr read_rds write_rds
#' @export
bs_fetch_results <- function(out_dir) {
  method_names <- list.dirs(out_dir, full.names = FALSE, recursive = FALSE) %>% discard(~ . == "")

  map_df(method_names, function(method_name) {
    method_folder <- paste0(out_dir, method_name)
    output_metrics_file <- paste0(method_folder, "/output_metrics.rds")
    output_models_file <- paste0(method_folder, "/output_models.rds")
    qsubhandle_file <- paste0(method_folder, "/qsubhandle.rds")

    if (!file.exists(output_metrics_file) && file.exists(qsubhandle_file)) {
      cat(method_name, ": Attempting to retrieve output from cluster: ", sep = "")
      data <- readr::read_rds(qsubhandle_file)
      grid <- data$grid
      qsub_handle <- data$qsub_handle
      num_tasks <- qsub_handle$num_tasks

      output <- PRISM::qsub_retrieve(
        qsub_handle,
        wait = FALSE
      )

      # suppressing inevitable warnings:
      # when the job is still partly running, calling qacct will generate a warning
      # when the job is finished, calling qstat will generate a warning
      suppressWarnings({
        qacct_out <- PRISM::qacct(qsub_handle)
        qstat_out <- PRISM::qstat_j(qsub_handle)
      })

      if (!is.null(output)) {
        cat("Output found! Saving output.\n", sep = "")

        output_succeeded <- TRUE
        outputs <- map_df(seq_along(output), function(grid_i) {
          out <- output[[grid_i]]
          if(length(out) != 1 || !is.na(out)) {
            out
          } else {
            qsub_error <- attr(out, "qsub_error")

            qsub_memory <- qsub_handle$memory %>% str_replace("G$", "") %>% as.numeric
            qacct_memory <- qacct_out$maxvmem %>% str_replace("GB$", "") %>% as.numeric

            if (!is.na(qacct_memory) && length(qacct_memory) > 0 && qacct_memory > qsub_memory) {
              qsub_error <- "Memory limit exceeded"
            }

            # mimic eval_ind format
            data_frame(
              method_name = data$method$name,
              method_short_name = data$method$short_name,
              error_message = qsub_error,
              repeat_i = grid$repeat_i[[grid_i]],
              fold_i = grid$fold_i[[grid_i]],
              group_sel = grid$group_sel[[grid_i]],
              grid_i
            )
          }
        })
      } else {
        error_message <-
          if (is.null(qstat_out) || nrow(qstat_out) > 0) {
            "job is still running"
          } else {
            "qsub_retrieve of results failed -- no output was produced, but job is not running any more"
          }

        cat("Output not found. ", error_message, ".\n", sep = "")
        output_succeeded <- FALSE

        NULL
      }

      if (output_succeeded) {
        if ("model" %in% colnames(outputs)) {
          models <- outputs$model
          model_ids <- map_chr(models, ~.$id)
          outputs <- outputs %>% slice(-model) %>% mutate(model_i = seq_len(n()), model_id = model_ids)
          readr::write_rds(models, output_models_file)
        }

        readr::write_rds(outputs, output_metrics_file)
      }

      gc()

      NULL

    } else {
      if (file.exists(output_metrics_file)) {
        cat(method_name, ": Output already present.\n", sep = "")
      } else {
        cat(method_name, ": No qsub file was found.\n", sep = "")
      }

      NULL
    }

  })
}



