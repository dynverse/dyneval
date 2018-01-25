#' A benchmark suite with which to run all the methods on the different tasks
#'
#' @param tasks A tibble of tasks.
#' @param task_group A grouping vector for the different tasks.
#' @param task_fold A fold index vector for the different tasks.
#' @param out_dir A folder in which to output intermediate and final results.
#' @param remote_dir A folder in which to store intermediate results in a remote directory when using the PRISM package.
#' @param timeout The number of seconds 1 method has to solve each of the tasks before a timeout is generated.
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
#' @param do_it_local Whether or not to run the benchmark suite locally (not recommended)
#' @param output_model Whether or not to return the outputted models
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
#'
#' @export
benchmark_suite_submit <- function(
  tasks,
  task_group,
  task_fold,
  out_dir,
  remote_dir,
  timeout = 120,
  methods = get_descriptions(as_tibble = TRUE),
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
  do_it_local = FALSE,
  output_model = FALSE
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
  mlr::configureMlr(show.learner.output = FALSE)
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

  ## save tasks to local file
  if (!do_it_local) {
    tasks_local_file <- paste0(out_dir, "/tasks_benchmark.rds")
    tasks_remote_file <- paste0(remote_dir, "/tasks_benchmark.rds")

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

    cat("Saving tasks file\n")
    PRISM:::mkdir_remote(path = out_dir, remote = "")
    readr::write_rds(tasks, tasks_local_file)

    cat("Moving tasks file to remote\n")
    PRISM:::mkdir_remote(path = remote_dir, remote = qsub_config$remote)
    PRISM:::rsync_remote("", tasks_local_file, qsub_config$remote, remote_dir)
  }

  ## run MBO for each method separately
  lapply (seq_len(nrow(methods)), function(methodi) {
    method <- dynutils::extract_row_to_list(methods, methodi)

    # determine where to store certain outputs
    method_folder <- paste0(out_dir, method$short_name)
    output_file <- paste0(method_folder, "/output.rds")
    qsubhandle_file <- paste0(method_folder, "/qsubhandle.rds")

    dir.create(method_folder, recursive = TRUE, showWarnings = FALSE)

    ## If no output or qsub handle exists yet
    if ((!file.exists(output_file) && !file.exists(qsubhandle_file)) || do_it_local) {
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

      # the function to run on the cluster
      qsub_x <- seq_len(nrow(grid))

      if (!do_it_local) {
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
          "method", "obj_fun", "design",
          "task_group", "task_fold", "tasks_remote_file",
          "num_cores", "metrics", "extra_metrics",
          "control", "grid",
          "learner", "timeout", "output_model",
          "num_folds"
        )

        # submit to the cluster
        qsub_handle <- PRISM::qsub_lapply(
          X = qsub_x,
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
          control, control, grid, qsub_handle
        )
        readr::write_rds(out, qsubhandle_file)

        invisible()
      } else {
        # run locally
        out <- pbapply::pblapply(
          X = qsub_x,
          FUN = benchmark_qsub_fun
        )
        out
      }
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
benchmark_qsub_fun <- function(grid_i) {
  fold_i <- grid[grid_i,]$fold_i
  group_sel <- grid[grid_i,]$group_sel
  repeat_i <- grid[grid_i,]$repeat_i

  if (!"tasks" %in% ls()) {
    tasks <- readr::read_rds(tasks_remote_file)
  }

  # create an objective function
  obj_fun <- make_obj_fun(method = method, metrics = metrics, extra_metrics = extra_metrics)

  ## start parameter optimisation
  parallelMap::parallelStartMulticore(cpus = num_cores, show.info = TRUE)

  if (num_folds != 1) {
    tune_train <- mlrMBO::mbo(
      obj_fun,
      learner = learner,
      design = design,
      control = control,
      show.info = TRUE,
      more.args = list(
        tasks = tasks[task_group == group_sel & task_fold != fold_i,],
        timeout = timeout,
        output_model = output_model
      )
    )
    design_test <- tune_train$opt.path$env$path %>% select(-one_of(metrics))
  } else {
    tune_train <- NULL
    design_test <- design
  }

  tune_test <- no_train_mlrMBO(
    obj_fun,
    learner = learner,
    design = design_test,
    control = control,
    show.info = TRUE,
    more.args = list(
      tasks = tasks[task_group == group_sel & task_fold == fold_i,],
      timeout = timeout,
      output_model = output_model
    )
  )
  parallelMap::parallelStop()

  list(tune_train = tune_train, tune_test = tune_test)
}

# Helper function for running just the initialisation of mlrMBO
no_train_mlrMBO <- function(fun, design, learner, control, show.info, more.args) {
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
#'
#' @importFrom PRISM qsub_retrieve qacct qstat_j
#' @importFrom readr read_rds write_rds
#' @export
benchmark_suite_retrieve <- function(out_dir) {
  method_names <- list.dirs(out_dir, full.names = FALSE, recursive = FALSE) %>% discard(~ . == "")
  map_df(method_names, function(method_name) {
    method_folder <- paste0(out_dir, method_name)
    output_file <- paste0(method_folder, "/output.rds")
    qsubhandle_file <- paste0(method_folder, "/qsubhandle.rds")

    if (file.exists(output_file)) {
      cat(method_name, ": Reading previous output\n", sep = "")
      readr::read_rds(output_file)
    } else if (!file.exists(output_file) && file.exists(qsubhandle_file)) {
      cat(method_name, ": Attempting to retrieve output from cluster: ", sep = "")
      data <- readr::read_rds(qsubhandle_file)
      qsub_handle <- data$qsub_handle
      num_tasks <- qsub_handle$num_tasks

      output <- PRISM::qsub_retrieve(
        qsub_handle,
        post_fun = function(rds_i, out_rds) {
          benchmark_suite_retrieve_helper(rds_i, out_rds, data)
        },
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
        outputs <- lapply(output, function(x) {
          if(length(x) == 1 && is.na(x)) {
            list(
              which_errored = list(TRUE),
              qsub_error = list(attr(x, "qsub_error"))
            )
          } else {
            # check to see whether all jobs failed
            all_errored <- all(x$individual_scores$error %>%  map_lgl(~ !is.null(.)))
            x$which_errored <- list(all_errored)
            x$qsub_error <- list(ifelse(all_errored, "all parameter settings errored", ""))
            x
          }
        }) %>% list_as_tibble()
      } else {
        error_message <-
          if (is.null(qstat_out) || nrow(qstat_out) > 0) {
            "job is still running"
          } else {
            "qsub_retrieve of results failed -- no output was produced, but job is not running any more"
          }

        cat("Output not found. ", error_message, ".\n", sep = "")
        output_succeeded <- FALSE
        outputs <- tibble(
          which_errored = list(rep(TRUE, num_tasks)),
          qsub_error = list(rep(error_message, num_tasks))
        )
      }

      out_rds <- outputs %>% mutate(
        task_id = seq_len(n()),
        method_name,
        qacct = map(seq_len(n()), ~ extract_row_to_list(qacct_out, .)),
        qstat = map(seq_len(n()), ~ qstat_out),
        qsub_handle = map(seq_len(n()), ~ qsub_handle)
      ) %>%
        select(method_name, task_id, which_errored, qsub_error, everything())

      if (output_succeeded) {
        readr::write_rds(out_rds, output_file)
      }

      out_rds
    } else {
      stop("Could not find an output.rds or qsubhandle.rds file in out_dir = ", sQuote(out_dir))
    }
  })
}

benchmark_suite_retrieve_helper <- function(rds_i, out_rds, data) {
  grid <- data$grid
  tasks <- data$tasks
  task_group <- data$task_group
  task_fold <- data$task_fold

  grid_i <- rds_i
  fold_i <- grid$fold_i[[rds_i]]
  group_sel <- grid$group_sel[[rds_i]]
  repeat_i <- grid$repeat_i[[rds_i]]
  design <- out_rds$design
  train_out <- out_rds$tune_train
  test_out <- out_rds$tune_test

  best_ind <- train_out$best.ind
  x_names <- names(train_out$x)
  y_names <- train_out$opt.path$y.names
  best <- list(
    param_index = best_ind,
    params = train_out$x,
    y_names = y_names,
    train_score = train_out$y,
    test_score = test_out$opt.path$env$path[best_ind,y_names]
  ) %>%
    list() %>%
    list_as_tibble() %>%
    mutate(
      grid_i, repeat_i, fold_i, group_sel
    )

  ## construct the global summary, 2 rows per parameter
  path_summary <- bind_rows(
    bind_cols(
      data_frame(type = "train", grid_i, repeat_i, fold_i, group_sel,
                 iteration_i = train_out$opt.path$env$dob,
                 time_total = train_out$opt.path$env$exec.time,
                 param_i = seq_along(time_total)),
      train_out$opt.path$env$path
    ),
    bind_cols(
      data_frame(type = "test", grid_i, repeat_i, fold_i, group_sel,
                 iteration_i = train_out$opt.path$env$dob,
                 time_total = test_out$opt.path$env$exec.time,
                 param_i = seq_along(time_total)),
      test_out$opt.path$env$path
    )
  ) %>% as.tibble
  par_set <- train_out$opt.path$par.set
  # par_set <- data$method$par_set

  for (parname in names(par_set$pars)) {
    parlen <- par_set$pars[[parname]]$len
    if (parlen > 1) {
      parlist <-
        path_summary %>%
        select(starts_with(parname)) %>%
        t %>%
        as.data.frame %>%
        as.list %>%
        set_names(NULL)
      path_summary <- path_summary %>%
        select(-starts_with(parname)) %>%
        mutate(!!parname := parlist)
    }
  }

  ## collect the scores per task individually
  individual_scores <- map_df(seq_along(train_out$opt.path$env$dob), function(param_i) {
    train_summary <- train_out$opt.path$env$extra[[param_i]]$.summary
    if (!is.null(train_summary)) {
      train_summary <- train_summary %>% mutate(fold_type = "train")
    }

    test_summary <- test_out$opt.path$env$extra[[param_i]]$.summary
    if (!is.null(test_summary)) {
      test_summary <- test_summary %>% mutate(fold_type = "test")
    }

    bind_rows(train_summary, test_summary) %>%
      mutate(grid_i, repeat_i, fold_i, group_sel, param_i,
             iteration_i = train_out$opt.path$env$dob[param_i])
  })

  lst(best, path_summary, individual_scores)
}

