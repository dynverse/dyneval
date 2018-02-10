#' A benchmark suite with which to run all the methods on the different tasks
#'
#' @param tasks A tibble of tasks.
#' @param task_group A grouping vector for the different tasks.
#' @param task_fold A fold index vector for the different tasks.
#' @param out_dir A folder in which to output intermediate and final results.
#' @param remote_dir A folder in which to store intermediate results in a remote directory when using the PRISM package.
#' @param optimisation_timeout The number of seconds a method is allowed to optimise parameters in a fold before generating a timeout.
#' @param methods A tibble of TI methods.
#' @param designs A names list of given designs data frames. Names must be equal to method short names.
#' @param metrics Which metrics to use;
#'   see \code{\link{calculate_metrics}} for a list of which metrics are available.
#' @param extra_metrics Extra metrics to calculate but not evaluate with.
#' @param num_cores The number of cores to allocate per mlr run.
#' @param memory The memory to allocate per core.
#' @param max_wall_time The maximum amount of time each fold is allowed to run.
#' @param execute_before Shell commands to execute before running R
#' @param r_module Which R module to use in gridegine
#' @param num_iterations The number of iterations to run.
#' @param num_init_params The number of initial parameters to evaluate.
#' @param num_repeats The number of times to repeat the mlr process, for each group and each fold.
#' @param output_model Whether or not to return the outputted models
#' @param verbose Whether or not to print extra information during parameter training
#'
#' @importFrom testthat expect_equal expect_is
#' @importFrom PRISM qsub_lapply override_qsub_config
#' @importFrom mlrMBO makeMBOControl setMBOControlTermination setMBOControlInfill makeMBOInfillCritDIB makeMBOInfillCritCB
#' @importFrom mlr makeLearner configureMlr
#' @importFrom ParamHelpers generateDesignOfDefaults generateDesign
#' @importFrom parallelMap parallelStartMulticore parallelStop
#' @importFrom randomForest randomForest
#' @importFrom pbapply pblapply
#' @importFrom emoa emoa_control
#' @importFrom readr read_rds write_rds
#' @importFrom dynmethods get_descriptions
#'
#' @export
benchmark_suite_submit <- function(
  tasks,
  task_group,
  task_fold,
  out_dir,
  remote_dir,
  optimisation_timeout = 24 * 3600,
  methods = dynmethods::get_descriptions(as_tibble = TRUE),
  designs = NULL,
  metrics = "correlation",
  extra_metrics = NULL,
  num_cores = 4,
  memory = "20G",
  max_wall_time = NULL,
  execute_before = NULL,
  r_module = "R",
  num_iterations = 20,
  num_init_params = 100,
  num_repeats = 1,
  output_model = FALSE,
  verbose = FALSE
) {
  testthat::expect_is(tasks, "tbl")
  testthat::expect_equal(nrow(tasks), length(task_group))
  testthat::expect_equal(nrow(tasks), length(task_fold))
  testthat::expect_is(methods, "tbl")

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
  mlr::configureMlr(show.learner.output = FALSE, on.learner.warning = "warn")
  learner <- mlr::makeLearner(
    "regr.randomForest",
    se.method = "jackknife",
    predict.type = "se",
    keep.inbag = TRUE
  )

  ## create a grid for each of the folds, groups, and repeats.
  ## one task will be created per row in this grid
  fold_ixs <- sort(unique(task_fold))
  num_folds <- length(fold_ixs)
  grid <- expand.grid(
    fold_i = fold_ixs,
    group_sel = sort(unique(task_group)),
    repeat_i = seq_len(num_repeats),
    stringsAsFactors = FALSE
  )

  ## Prepare for remote execution; create a qsub config
  qsub_config <- PRISM::override_qsub_config(
    wait = FALSE,
    remove_tmp_folder = FALSE,
    stop_on_error = FALSE,
    verbose = FALSE,
    num_cores = num_cores,
    memory = memory,
    max_wall_time = max_wall_time,
    r_module = r_module,
    execute_before = execute_before,
    local_tmp_path = out_dir,
    remote_tmp_path = remote_dir
  )

  ## Save tasks to local file
  tasks_local_file <- paste0(out_dir, "/tasks_benchmark.rds")
  tasks_remote_file <- paste0(remote_dir, "/tasks_benchmark.rds")

  cat("Saving tasks file\n")
  PRISM:::mkdir_remote(path = out_dir, remote = "")
  readr::write_rds(tasks, tasks_local_file)

  cat("Moving tasks file to remote\n")
  PRISM:::mkdir_remote(path = remote_dir, remote = qsub_config$remote)
  PRISM:::rsync_remote("", tasks_local_file, qsub_config$remote, remote_dir)

  ## run MBO for each method separately
  lapply (seq_len(nrow(methods)), function(methodi) {
    method <- dynutils::extract_row_to_list(methods, methodi)

    # determine where to store certain outputs
    method_folder <- paste0(out_dir, "/", method$short_name)
    output_file <- paste0(method_folder, "/output.rds")
    qsubhandle_file <- paste0(method_folder, "/qsubhandle.rds")

    PRISM:::mkdir_remote(path = method_folder, remote = "")

    ## If no output or qsub handle exists yet
    if (!file.exists(output_file) && !file.exists(qsubhandle_file)) {
      cat("Submitting ", method$name, "\n", sep="")

      # generate initial parameters
      if (is.null(designs)) {
        design <- bind_rows(
          ParamHelpers::generateDesignOfDefaults(method$par_set),
          ParamHelpers::generateDesign(n = num_init_params, par.set = method$par_set)
        )
      } else {
        design <- designs[[method$short_name]]
      }

      # set parameters for the cluster
      qsub_config_method <- PRISM::override_qsub_config(
        qsub_config = qsub_config,
        name = paste0("D_", method$short_name),
        local_tmp_path = paste0(method_folder, "/r2gridengine")
      )

      # which packages to load on the cluster
      qsub_packages <- c("dplyr", "purrr", "dyneval", "mlrMBO", "parallelMap", "readr")

      # which data objects will need to be transferred to the cluster
      qsub_environment <-  c(
        "method", "design", "task_group", "task_fold", "tasks_remote_file",
        "num_cores", "metrics", "extra_metrics", "control", "grid",
        "learner", "optimisation_timeout", "output_model", "num_folds", "verbose"
      )

      # submit to the cluster
      qsub_handle <- PRISM::qsub_lapply(
        X = seq_len(nrow(grid)),
        object_envir = environment(),
        qsub_environment = qsub_environment,
        qsub_packages = qsub_packages,
        qsub_config = qsub_config_method,
        FUN = benchmark_qsub_fun
      )

      # save data and handle to RDS file
      out <- lst(
        method, design,
        task_group, task_fold,
        num_cores, metrics, extra_metrics,
        control, grid, qsub_handle,
        tasks_local_file, tasks_remote_file
      )
      readr::write_rds(out, qsubhandle_file)

      NULL
    }
  })
}

#' Helper function for benchmark suite
#'
#' @param grid_i Benchmark config index
#'
#' @importFrom readr read_rds
#' @importFrom parallelMap parallelStartMulticore parallelStop
#' @importFrom mlrMBO mbo
#' @importFrom mlr configureMlr
benchmark_qsub_fun <- function(grid_i) {
  fold_i <- grid[grid_i,]$fold_i
  group_sel <- grid[grid_i,]$group_sel
  repeat_i <- grid[grid_i,]$repeat_i

  working_dir <- getwd()

  # configure mlr
  mlr::configureMlr(show.learner.output = FALSE, on.learner.warning = "warn")

  if (!"tasks" %in% ls()) {
    tasks <- readr::read_rds(tasks_remote_file)
  }

  # create an objective function
  obj_fun <- make_obj_fun(method = method, metrics = metrics, extra_metrics = extra_metrics, verbose = verbose)

  if (num_folds != 1) {
    ## create a folder to save the intermediate mbo files in
    save_file_path <- paste0(working_dir, "/mlrmbo")
    dir.create(save_file_path, showWarnings = FALSE, recursive = TRUE)

    ## configure intermediate output
    control_train <- control
    control_train$save.file.path <- paste0(save_file_path, "/mlr_progress_", grid_i, ".RData")
    control_train$save.on.disk.at <- seq(0, control_train$iters+1, by = 1)

    ## start parallellisation
    parallelMap::parallelStartMulticore(cpus = num_cores, show.info = TRUE)

    ## Start training. If a optimisation_timeout is reached, the training will stop prematurely.
    train_out_poss_timeout <- dynutils::eval_with_timeout(
      timeout = optimisation_timeout,
      expr = {
        mlrMBO::mbo(
          obj_fun,
          learner = learner,
          design = design,
          control = control_train,
          show.info = TRUE,
          more.args = list(
            tasks = tasks[task_group == group_sel & task_fold != fold_i,],
            output_model = output_model
          )
        )
      }
    )

    ## stop parallellisation
    parallelMap::parallelStop()

    ## Read the last saved state
    train_out <- mlrMBO::mboFinalize(control_train$save.file.path)

    ## Extract all parameters from the training, to be run on the test datasets
    design_test <- train_out$opt.path$env$path %>% select(-one_of(metrics))
  } else {
    train_out <- NULL
    design_test <- design
  }

  if (num_folds != 1) {
    ## start parallellisation
    parallelMap::parallelStartMulticore(cpus = num_cores, show.info = TRUE)
    mc_cores <- 1
  } else {
    mc_cores = num_cores
  }

  test_out <- no_train_mlrmbo(
    obj_fun,
    learner = learner,
    design = design_test,
    control = control,
    show.info = TRUE,
    more.args = list(
      tasks = tasks[task_group == group_sel & task_fold == fold_i,],
      output_model = output_model,
      mc_cores = num_cores
    )
  )

  if (num_folds != 1) {
    ## stop parallellisation
    parallelMap::parallelStop()
  }

  # get iterations
  if (!is.null(train_out)) {
    iterations <- train_out$opt.path$env$dob
  } else {
    iterations <- test_out$opt.path$env$dob
  }

  # get each individual evaluation (model, params, performance, timings) and put it in a tibble
  map_df(seq_len(nrow(design_test)), function(param_i) {
    bind_rows(
      summarise_mlr_output(fold_type = "train", train_out, param_i, design_test, output_model, grid, grid_i, iterations),
      summarise_mlr_output(fold_type = "test", test_out, param_i, design_test, output_model, grid, grid_i, iterations)
    )
  })
}

summarise_mlr_output <- function(fold_type, mlr_out, param_i, design, output_model, grid, grid_i, iterations) {
  if (!is.null(mlr_out)) {
    summary <- mlr_out$opt.path$env$extra[[param_i]]$.summary

    if (!is.null(summary)) {
      summary <- summary %>%
        mutate(
          fold_type,
          param_row = list(design[param_i,]),
          grid_i,
          repeat_i = grid$repeat_i[grid_i],
          fold_i = grid$fold_i[grid_i],
          group_sel = grid$group_sel[grid_i],
          param_i,
          iteration_i = iterations[param_i],
          error_message = sapply(error, function(err) ifelse(is.null(err), "", err$message))
        ) %>%
        select(-error)

      # save models, if asked for
      if (output_model) {
        summary$model <- mlr_out$opt.path$env$extra[[param_i]]$.models
      }
    }

    summary
  } else {
    NULL
  }
}

# Helper function for running just the initialisation of mlrMBO
no_train_mlrmbo <- function(fun, design, learner, control, show.info, more.args) {
  requireNamespace("mlrMBO")
  opt.problem <- mlrMBO:::initOptProblem(
    fun = fun,
    design = design,
    learner = learner,
    control = control,
    show.info = show.info,
    more.args = more.args
  )
  opt.state <- mlrMBO:::makeOptState(opt.problem)
  mlrMBO:::evalMBODesign.OptState(opt.state)
  mlrMBO:::finalizeMboLoop(opt.state)
  opt.result <- mlrMBO:::getOptStateOptResult(opt.state)
  mbo.result <- mlrMBO:::makeMBOResult.OptState(opt.state)
  mbo.result
}


#' Downloading and processing the results of the benchmark jobs
#'
#' @param out_dir The folder in which to output intermediate and final results.
#' @param return_outputs Whether or not to return the outputs, or only retrieve them from the cluster.
#'
#' @importFrom PRISM qsub_retrieve qacct qstat_j
#' @importFrom readr read_rds write_rds
#' @export
benchmark_suite_retrieve <- function(out_dir, return_outputs = TRUE) {
  method_names <- list.dirs(out_dir, full.names = FALSE, recursive = FALSE) %>% discard(~ . == "")
  map_df(method_names, function(method_name) {
    method_folder <- paste0(out_dir, method_name)
    output_file <- paste0(method_folder, "/output.rds")
    qsubhandle_file <- paste0

    if (return_outputs) {
      if (file.exists(output_file)) {
        cat(method_name, ": Reading previous output\n", sep = "")
        readr::read_rds(output_file)
      } else {
        cat(method_name, ": Output not found, skipping\n", sep = "")
        NULL
      }
    } else
    {
      if (!file.exists(output_file) && file.exists(qsubhandle_file)) {
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
          readr::write_rds(outputs, output_file)
        }

        gc()

        NULL

      } else {
        if (file.exists(output_file)) {
          cat(method_name, ": Output already present", sep = "")
        } else {
          cat(method_name, ": No qsub file was found", sep = "")
        }

        NULL
      }

    }
  })
}



