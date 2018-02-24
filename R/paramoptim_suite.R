#' A parameter optimisation suite
#'
#' @param task_ids The ids of the tasks to be used in the evaluation.
#' @param local_tasks_folder A local folder containing all of the tasks to be evaluated on.
#' @param remote_tasks_folder A remote folder containing all of the tasks to be evaluated on.
#' @param methods A tibble of TI methods to use, see \code{\link[dynmethods]{get_descriptions}}.
#' @param timeout_paramoptim The parameter optimisation timeout
#' @param max_memory_per_core The maximum amount of memory each core is allowed to use
#' @param num_cores The number of cores to use.
#' @param metrics Which metrics to evaluate; see \code{\link{calculate_metrics}} for a list of which metrics are available.
#' @param num_repeats The number of times to repeat the evaluation.
#' @param num_iterations The number of iterations to run.
#' @param num_init_params The number of initial parameters to evaluate.
#' @param local_output_folder A folder in which to output intermediate and final results.
#' @param remote_output_folder A folder in which to store intermediate results in a remote directory when using the PRISM package.
#' @param execute_before Shell commands to execute before running R.
#' @param verbose Whether or not to print extra information.
#'
#' @importFrom PRISM qsub_lapply override_qsub_config
#' @importFrom pbapply pblapply
#' @importFrom readr read_rds write_rds
#' @importFrom testthat expect_equal expect_is
#' @importFrom mlrMBO makeMBOControl setMBOControlTermination setMBOControlInfill makeMBOInfillCritDIB makeMBOInfillCritCB
#' @importFrom mlr makeLearner configureMlr
#' @importFrom ParamHelpers generateDesignOfDefaults generateDesign
#' @importFrom parallelMap parallelStartMulticore parallelStop
#' @importFrom randomForest randomForest
#' @importFrom emoa emoa_control
#' @importFrom dynmethods get_descriptions
#'
#' @export
paramoptim_submit <- function(
  task_ids,
  local_tasks_folder,
  remote_tasks_folder,
  methods,
  timeout_paramoptim,
  max_memory_per_core,
  num_cores,
  metrics,
  num_repeats,
  num_iterations,
  num_init_params,
  local_output_folder,
  remote_output_folder,
  execute_before = NULL,
  verbose = FALSE
) {
  paramoptim_submit_check(
    local_tasks_folder,
    remote_tasks_folder,
    task_ids,
    methods,
    timeout_paramoptim,
    max_memory_per_core,
    metrics,
    num_repeats,
    num_iterations,
    num_init_params,
    local_output_folder,
    remote_output_folder,
    execute_before,
    verbose
  )

  ##################### MBO CONFIG ################################
  ## set settings for MBO parameter optimisation
  control <- mlrMBO::makeMBOControl(
    n.objectives = length(metrics),
    propose.points = num_cores,
    y.name = metrics
  ) %>%
    mlrMBO::setMBOControlTermination(
      iters = num_iterations
    )

  # Set infill criterion
  if (length(metrics) == 1) {
    control <- control %>%
      mlrMBO::setMBOControlInfill(mlrMBO::makeMBOInfillCritCB())
  } else {
    control <- control %>%
      mlrMBO::setMBOControlInfill(mlrMBO::makeMBOInfillCritDIB())
  }

  ## create learner for predicting performance of new params
  learner <- mlr::makeLearner(
    "regr.randomForest",
    se.method = "jackknife",
    predict.type = "se",
    keep.inbag = TRUE
  )

  ##################### PRISM CONFIG ################################
  ## prepare for remote execution; create a qsub config
  qsub_config <- PRISM::override_qsub_config(
    wait = FALSE,
    remove_tmp_folder = FALSE,
    stop_on_error = FALSE,
    verbose = FALSE,
    num_cores = num_cores,
    memory = max_memory_per_core,
    max_wall_time = timeout_paramoptim,
    r_module = NULL,
    execute_before = execute_before,
    local_tmp_path = local_output_folder,
    remote_tmp_path = remote_output_folder
  )

  ## run evaluation for each method separately
  map(seq_len(nrow(methods)), function(methodi) {
    method <- dynutils::extract_row_to_list(methods, methodi)

    # determine where to store certain outputs
    method_folder <- paste0(local_output_folder, "/", method$short_name)
    output_file <- paste0(method_folder, "/output.rds")
    qsubhandle_file <- paste0(method_folder, "/qsubhandle.rds")

    PRISM:::mkdir_remote(path = method_folder, remote = "")

    ## If no output or qsub handle exists yet
    if (!file.exists(output_file) && !file.exists(qsubhandle_file)) {
      cat("Submitting ", method$name, "\n", sep="")

      # fetch parameters
      parms <- parameters[[method$short_name]]

      ## create a grid for each of the tasks, paramsets, and repeats
      grid <- crossing(
        repeat_i = seq_len(num_repeats)
      )

      # generate initial parameters
      design <- bind_rows(
        ParamHelpers::generateDesignOfDefaults(method$par_set),
        ParamHelpers::generateDesign(n = num_init_params, par.set = method$par_set)
      )

      # set parameters for the cluster
      qsub_config_method <-
        PRISM::override_qsub_config(
          qsub_config = qsub_config,
          name = paste0("D_", method$short_name),
          local_tmp_path = paste0(method_folder, "/r2gridengine")
        )

      # which packages to load on the cluster
      qsub_packages <- c("dplyr", "purrr", "dyneval", "mlrMBO", "parallelMap", "readr")

      # which data objects will need to be transferred to the cluster
      qsub_environment <-  c(
        "grid", "remote_tasks_folder", "remote_output_folder", "parms", "method", "metrics", "verbose",
        "num_cores", "control", "learner", "design", "task_ids"
      )

      # submit to the cluster
      qsub_handle <- PRISM::qsub_lapply(
        X = seq_len(nrow(grid)),
        object_envir = environment(),
        qsub_environment = qsub_environment,
        qsub_packages = qsub_packages,
        qsub_config = qsub_config_method,
        FUN = paramoptim_qsub_fun
      )

      # save data and handle to RDS file
      metadata <- lst(
        local_tasks_folder,
        remote_tasks_folder,
        local_output_folder,
        remote_output_folder,
        task_ids,
        method,
        timeout_paramoptim,
        max_memory_per_core,
        num_cores,
        metrics,
        num_repeats,
        num_iterations,
        num_init_params,
        grid,
        control,
        learner,
        method_folder,
        output_file,
        qsubhandle_file,
        parms,
        qsub_handle
      )
      readr::write_rds(metadata, qsubhandle_file)

      NULL
    }
  })

  invisible()
}

#' @importFrom testthat expect_equal expect_is
#' @importFrom ParamHelpers dfRowToList
paramoptim_submit_check <- function(
  local_tasks_folder,
  remote_tasks_folder,
  task_ids,
  methods,
  timeout_paramoptim,
  max_memory_per_core,
  metrics,
  num_repeats,
  num_iterations,
  num_init_params,
  local_output_folder,
  remote_output_folder,
  execute_before,
  verbose
) {
  # check tasks
  # TODO: check whether all tasks are present, local and remote

  # check methods
  testthat::expect_is(methods, "tbl")

  # check timeout_per_execution
  testthat::expect_is(timeout_per_execution, "numeric")

  # check max_memory_per_execution
  testthat::expect_is(max_memory_per_execution, "character")
  testthat::expect_match(max_memory_per_execution, "[0-9]+G")
}


#' Helper function for paramoptim suite
#'
#' @param grid_i paramoptim config index
paramoptim_qsub_fun <- function(grid_i) {
  # call helper function
  paramoptim_run_evaluation(
    grid,
    grid_i,
    task_ids,
    control,
    learner,
    design,
    remote_tasks_folder,
    parms,
    method,
    metrics,
    num_cores,
    verbose
  )
}

#' @importFrom readr read_rds
#' @importFrom ParamHelpers dfRowToList trafoValue
#' @importFrom parallelMap parallelStartMulticore parallelStop
#' @importFrom mlrMBO mbo
#' @importFrom mlr configureMlr
paramoptim_run_evaluation <- function(
  grid,
  grid_i,
  task_ids,
  control,
  learner,
  design,
  remote_tasks_folder,
  parms,
  method,
  metrics,
  num_cores,
  verbose
) {
  # configure mlr
  mlr::configureMlr(show.learner.output = FALSE, on.learner.warning = "warn")

  # read tasks
  tasks <- map_df(paste0(remote_tasks_folder, "/", task_ids, ".rds"), readr::read_rds)

  # create an objective function
  obj_fun <- make_obj_fun(method = method, metrics = metrics, extra_metrics = NULL, verbose = verbose)

  ## create a folder to save the intermediate mbo files in
  save_file_path <- paste0(working_dir, "/mlrmbo")
  dir.create(save_file_path, showWarnings = FALSE, recursive = TRUE)

  ## configure intermediate output
  control_train <- control
  control_train$save.file.path <- paste0(save_file_path, "/mlr_progress_", grid_i, ".RData")
  control_train$save.on.disk.at <- seq(0, control_train$iters+1, by = 1)

  ## start parallellisation
  parallelMap::parallelStartMulticore(cpus = num_cores, show.info = TRUE)

  ## Start training, and write intermediate results frequently to file
  mlrMBO::mbo(
    obj_fun,
    learner = learner,
    design = design,
    control = control_train,
    show.info = TRUE,
    more.args = list(
      tasks = tasks,
      output_model = FALSE
    )
  )

  ## stop parallellisation
  parallelMap::parallelStop()

  invisible()
}


#' Fetch the results of the paramoptim jobs from the cluster.
#'
#' @param local_output_folder The folder in which to output intermediate and final results.
#'
#' @importFrom PRISM qsub_retrieve qacct qstat_j
#' @importFrom readr read_rds write_rds
#' @export
paramoptim_fetch_results <- function(local_output_folder) {
  method_names <- list.dirs(local_output_folder, full.names = FALSE, recursive = FALSE) %>% discard(~ . == "")

  # process each method separately
  map(method_names, function(method_name) {
    # method_folder <- paste0(local_output_folder, "/", method_name)
    # output_metrics_file <- paste0(method_folder, "/output_metrics.rds")
    # output_models_file <- paste0(method_folder, "/output_models.rds")
    # qsubhandle_file <- paste0(method_folder, "/qsubhandle.rds")
    #
    # # if the output has not been processed yet, but a qsub handle exists,
    # # attempt to fetch the results from the cluster
    # if (!file.exists(output_metrics_file) && file.exists(qsubhandle_file)) {
    #
    #   cat(method_name, ": Attempting to retrieve output from cluster: ", sep = "")
    #   metadata <- readr::read_rds(qsubhandle_file)
    #   grid <- metadata$grid
    #   qsub_handle <- metadata$qsub_handle
    #   num_tasks <- qsub_handle$num_tasks
    #
    #   # attempt to retrieve results; return NULL if job is still busy or has failed
    #   output <- PRISM::qsub_retrieve(
    #     qsub_handle,
    #     wait = FALSE
    #   )
    #
    #
    #
    #
    #   # ## Read the last saved state
    #   # train_out <- mlrMBO::mboFinalize(control_train$save.file.path)
    #   # if (!is.null(mlr_out)) {
    #   #   summary <- mlr_out$opt.path$env$extra[[param_i]]$.summary
    #   #
    #   #   if (!is.null(summary)) {
    #   #     summary <- summary %>%
    #   #       mutate(
    #   #         fold_type,
    #   #         param_row = list(design[param_i,]),
    #   #         grid_i,
    #   #         repeat_i = grid$repeat_i[grid_i],
    #   #         param_i,
    #   #         iteration_i = iterations[param_i],
    #   #         error_message = sapply(error, function(err) ifelse(is.null(err), "", err$message))
    #   #       ) %>%
    #   #       select(-error)
    #   #   }
    #   #
    #   #   summary
    #   # } else {
    #   #   NULL
    #   # }
    #
    #   if (!is.null(output)) {
    #     cat("Output found! Saving output.\n", sep = "")
    #
    #     suppressWarnings({
    #       qacct_out <- PRISM::qacct(qsub_handle)
    #     })
    #
    #     # process each job separately
    #     outputs <- map_df(seq_len(nrow(grid)), function(grid_i) {
    #       out <- output[[grid_i]]
    #
    #       if (length(out) != 1 || !is.na(out)) {
    #         # hooray, the paramoptim suite ran fine!
    #         out
    #
    #       } else {
    #         # if qacct is empty or the correct taskid cannot be found, then
    #         # this job never ran
    #         if (is.null(qacct_out) || !any(qacct_out$taskid == grid_i)) {
    #           qsub_error <- "Job cancelled by user"
    #         } else {
    #           qacct_filt <- qacct_out %>% filter(taskid == grid_i) %>% arrange(desc(row_number_i)) %>% slice(1)
    #
    #           qacct_memory <- qacct_filt$maxvmem %>% str_replace("GB$", "") %>% as.numeric
    #           qacct_exit_status <- qacct_filt$exit_status %>% str_replace(" +", " ")
    #           qacct_exit_status_number <- qacct_exit_status %>% str_replace(" .*", "")
    #           qacct_user_time <- qacct_filt$ru_stime %>% str_replace("s$", "") %>% as.numeric
    #
    #           qsub_memory <- qsub_handle$memory %>% str_replace("G$", "") %>% as.numeric
    #           qsub_user_time <- qsub_handle$max_wall_time
    #
    #           qsub_error <-
    #             #if (qacct_memory > qsub_memory) {
    #             if (qacct_exit_status_number %in% c("134", "139")) {
    #               "Memory limit exceeded"
    #               #} else if (qacct_user_time > qacct_user_time) {
    #             } else if (qacct_exit_status_number %in% c("137", "140", "9", "64")) {
    #               "Time limit exceeded"
    #             } else if (qacct_exit_status_number != "0") {
    #               qacct_exit_status
    #             } else {
    #               attr(out, "qsub_error")
    #             }
    #         }
    #
    #         # create a method that will just return the error generated by qsub
    #         method_failer <- metadata$method
    #         method_failer$run_fun <- function(...) {
    #           stop(qsub_error)
    #         }
    #         formals(method_failer$run_fun) <- formals(metadata$method$run_fun)
    #
    #         # "rerun" the evaluation, in order to generate the expected output except with
    #         # the default fail-scores for each of the metrics
    #         out <- paramoptim_run_evaluation(
    #           grid = grid,
    #           grid_i = grid_i,
    #           remote_tasks_folder = metadata$local_tasks_folder,
    #           parms = metadata$parms,
    #           method = method_failer,
    #           metrics = metadata$metrics,
    #           verbose = FALSE
    #         )
    #       }
    #     })
    #
    #     # save models separately
    #     models <- outputs$model
    #     model_ids <- map_chr(models, function(model) {
    #       if (!is.null(model)) {
    #         model$id
    #       } else {
    #         NA
    #       }
    #     })
    #     models <- models %>% setNames(model_ids)
    #     outputs <- outputs %>% select(-model) %>% mutate(model_i = seq_len(n()), model_id = model_ids)
    #     readr::write_rds(models, output_models_file)
    #
    #     # save output
    #     readr::write_rds(outputs, output_metrics_file)
    #
    #   } else {
    #     # the job is probably still running
    #     suppressWarnings({
    #       qstat_out <- PRISM::qstat_j(qsub_handle)
    #     })
    #
    #     error_message <-
    #       if (is.null(qstat_out) || nrow(qstat_out) > 0) {
    #         "job is still running"
    #       } else {
    #         "qsub_retrieve of results failed -- no output was produced, but job is not running any more"
    #       }
    #
    #     cat("Output not found. ", error_message, ".\n", sep = "")
    #   }
    #
    #   NULL
    # } else {
    #   if (file.exists(output_metrics_file)) {
    #     cat(method_name, ": Output already present.\n", sep = "")
    #   } else {
    #     cat(method_name, ": No qsub file was found.\n", sep = "")
    #   }
    #   NULL
    # }

  })

  # return nothing
  invisible()
}



#' Gather and bind the results of the paramoptim jobs
#'
#' @param local_output_folder The folder in which to output intermediate and final results.
#' @param load_models Whether or not to load the models as well.
#'
#' @importFrom readr read_rds
#' @export
paramoptim_bind_results <- function(local_output_folder, load_models = FALSE) {
  method_names <- list.dirs(local_output_folder, full.names = FALSE, recursive = FALSE) %>% discard(~ . == "")

  # process each method separately
  as_tibble(map_df(method_names, function(method_name) {
    # method_folder <- paste0(local_output_folder, method_name)
    # output_metrics_file <- paste0(method_folder, "/output_metrics.rds")
    # output_models_file <- paste0(method_folder, "/output_models.rds")
    #
    # if (file.exists(output_metrics_file)) {
    #   cat(method_name, ": Reading previous output\n", sep = "")
    #   output <- readr::read_rds(output_metrics_file)
    #
    #   # read models, if requested
    #   if (load_models && file.exists(output_models_file)) {
    #     models <- readr::read_rds(output_models_file)
    #     output$model <- models
    #   }
    #
    #   output
    # } else {
    #   cat(method_name, ": Output not found, skipping\n", sep = "")
    #   NULL
    # }
  }))
}

