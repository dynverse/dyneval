#' Fetch the results of the benchmark jobs from the cluster.
#'
#' @param out_dir The folder in which to output intermediate and final results.
#'
#' @importFrom PRISM qsub_retrieve qacct qstat_j
#' @importFrom readr read_rds write_rds
#' @export
bs_fetch_results <- function(out_dir) {
  method_names <- list.dirs(out_dir, full.names = FALSE, recursive = FALSE) %>% discard(~ . == "")

  # process each method separately
  map(method_names, function(method_name) {
    method_folder <- paste0(out_dir, method_name)
    output_metrics_file <- paste0(method_folder, "/output_metrics.rds")
    output_models_file <- paste0(method_folder, "/output_models.rds")
    qsubhandle_file <- paste0(method_folder, "/qsubhandle.rds")

    # if the output has not been processed yet, but a qsub handle exists,
    # attempt to fetch the results from the cluster
    if (!file.exists(output_metrics_file) && file.exists(qsubhandle_file)) {

      cat(method_name, ": Attempting to retrieve output from cluster: ", sep = "")
      data <- readr::read_rds(qsubhandle_file)
      grid <- data$grid
      qsub_handle <- data$qsub_handle
      num_tasks <- qsub_handle$num_tasks

      # attempt to retrieve results; return NULL if job is still busy or has failed
      output <- PRISM::qsub_retrieve(
        qsub_handle,
        wait = FALSE
      )

      if (!is.null(output)) {
        cat("Output found! Saving output.\n", sep = "")

        # process each job separately
        outputs <- map_df(seq_along(output), function(grid_i) {
          out <- output[[grid_i]]

          if (length(out) != 1 || !is.na(out)) {
            # hooray, the benchmark suite ran fine!
            out

          } else {
            # a subtask has errored, will mimic regular output
            qsub_error <- attr(out, "qsub_error")

            suppressWarnings({
              qacct_out <- PRISM::qacct(qsub_handle)
            })

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

        if ("model" %in% colnames(outputs)) {
          models <- outputs$model
          model_ids <- map_chr(models, function(model) {
            if (!is.null(model)) {
              model$id
            } else {
              NA
            }
          })
          models <- models %>% setNames(model_ids)
          outputs <- outputs %>% select(-model) %>% mutate(model_i = seq_len(n()), model_id = model_ids)
          readr::write_rds(models, output_models_file)
        }

        # save output
        readr::write_rds(outputs, output_metrics_file)

      } else {
        # the whole job has errored, will mimic regular output

        suppressWarnings({
          qstat_out <- PRISM::qstat_j(qsub_handle)
        })

        error_message <-
          if (is.null(qstat_out) || nrow(qstat_out) > 0) {
            "job is still running"
          } else {
            "qsub_retrieve of results failed -- no output was produced, but job is not running any more"
          }

        cat("Output not found. ", error_message, ".\n", sep = "")
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

  # return nothing
  invisible()
}



