#' A benchmark suite with which to run all the methods on the different tasks
#'
#' @param tasks A tibble of tasks.
#' @param task_group A grouping vector for the different tasks.
#' @param task_fold A fold index vector for the different tasks.
#' @param out_dir The folder in which to output intermediate and final results.
#' @param timeout The number of seconds 1 method has to solve each of the tasks before a timeout is generated.
#' @param methods A tibble of TI methods.
#' @param metrics Which metrics to use;
#'   see \code{\link{calculate_metrics}} for a list of which metrics are available.
#' @param num_cores The number of cores to allocate per mlr run.
#' @param memory The memory to allocate per core.
#' @param max_wall_time The maximum amount of time each fold is allowed to run.
#' @param num_iterations The number of iterations to run.
#' @param num_init_params The number of initial parameters to evaluate.
#' @param num_repeats The number of times to repeat the mlr process, for each group and each fold.
#' @param save_r2g_to_outdir Save the r2gridengine output to \code{out_dir} instead of the default \code{local_tmp_path}.
#' @param do_it_local Whether or not to run the benchmark suite locally (not recommended)
#'
#' @importFrom testthat expect_equal expect_is
#' @importFrom PRISM qsub_lapply override_qsub_config
#' @importFrom mlrMBO makeMBOControl setMBOControlTermination setMBOControlInfill makeMBOInfillCritDIB makeMBOInfillCritCB
#' @importFrom mlr makeLearner configureMlr
#' @importFrom ParamHelpers generateDesignOfDefaults generateDesign
#' @importFrom parallelMap parallelStartMulticore parallelStop
#' @importFrom randomForest randomForest
#' @importFrom pbapply pblapply
#'
#' @export
benchmark_suite_submit <- function(
  tasks,
  task_group,
  task_fold,
  out_dir,
  timeout = 120,
  methods = get_descriptions(as_tibble = TRUE),
  metrics = c("auc_R_nx", "robbie_network_score"),
  num_cores = 4,
  memory = "20G",
  max_wall_time = NULL,
  num_iterations = 20,
  num_init_params = 100,
  num_repeats = 1,
  save_r2g_to_outdir = FALSE,
  do_it_local = FALSE
) {
  testthat::expect_is(tasks, "tbl")
  testthat::expect_equal(nrow(tasks), length(task_group))
  testthat::expect_equal(nrow(tasks), length(task_fold))
  testthat::expect_is(methods, "tbl")

  ## set settings for MBO parameter optimisation
  control_train <- mlrMBO::makeMBOControl(
    n.objectives = length(metrics),
    propose.points = num_cores,
    y.name = metrics
  ) %>%
    mlrMBO::setMBOControlTermination(
      iters = num_iterations
    )

  # Set infill criterion
  if (length(metrics) == 1) {
    control_train <- control_train %>%
      mlrMBO::setMBOControlInfill(mlrMBO::makeMBOInfillCritCB())
  } else {
    control_train <- control_train %>%
      mlrMBO::setMBOControlInfill(mlrMBO::makeMBOInfillCritDIB())
  }

  # construct control for test phase
  control_test <- control_train %>% mlrMBO::setMBOControlTermination(
    iters = 1
  )
  control_test$propose.points <- 1

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
  grid <- expand.grid(
    fold_i = sort(unique(task_fold)),
    group_sel = sort(unique(task_group)),
    repeat_i = seq_len(num_repeats),
    stringsAsFactors = FALSE
  )

  ## run MBO for each method separatelly
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

      # create an objective function
      obj_fun <- make_obj_fun(method = method, metrics = metrics)

      # generate initial parameters
      design <- bind_rows(
        ParamHelpers::generateDesignOfDefaults(method$par_set),
        ParamHelpers::generateDesign(n = num_init_params, par.set = method$par_set)
      )

      # the function to run on the cluster
      qsub_x <- seq_len(nrow(grid))

      qsub_fun <- function(grid_i) {
        fold_i <- grid[grid_i,]$fold_i
        group_sel <- grid[grid_i,]$group_sel
        repeat_i <- grid[grid_i,]$repeat_i

        ## start parameter optimisation
        # TODO: If the models should be outputted, change the output_model
        # and process the output

        parallelMap::parallelStartMulticore(cpus = num_cores, show.info = TRUE)
        tune_train <- mlrMBO::mbo(
          obj_fun,
          learner = learner,
          design = design,
          control = control_train,
          show.info = TRUE,
          more.args = list(
            tasks = tasks[task_group == group_sel & task_fold != fold_i,],
            timeout = timeout,
            output_model = FALSE #"models/"
          )
        )
        tune_test <- mlrMBO::mbo(
          obj_fun,
          learner = learner,
          design = tune_train$opt.path$env$path %>% select(-one_of(metrics)),
          control = control_test,
          show.info = TRUE,
          more.args = list(
            tasks = tasks[task_group == group_sel & task_fold == fold_i,],
            timeout = timeout,
            output_model = FALSE #"models/"
          )
        )
        parallelMap::parallelStop()

        list(design = design, tune_train = tune_train, tune_test = tune_test)
      }

      if (!do_it_local) {
        # set parameters for the cluster
        qsub_config <- PRISM::override_qsub_config(
          wait = FALSE,
          num_cores = num_cores,
          memory = memory,
          name = paste0("D_", method$short_name),
          remove_tmp_folder = FALSE,
          stop_on_error = FALSE,
          verbose = FALSE,
          max_wall_time = max_wall_time,
          execute_before = "export R_MAX_NUM_DLLS=300"
        )

        # if the cluster data needs to be saved to dyneval output folder
        if (save_r2g_to_outdir) {
          qsub_config$local_tmp_path <- paste0(method_folder, "/r2gridengine")
          dir.create(qsub_config$local_tmp_path, recursive = T, showWarnings = F)
        }

        # which packages to load on the cluster
        qsub_packages <- c("dplyr", "purrr", "dyneval", "mlrMBO", "parallelMap")

        # which data objects will need to be transferred to the cluster
        qsub_environment <-  c(
          "method", "obj_fun", "design",
          "tasks", "task_group", "task_fold",
          "num_cores", "metrics",
          "control_train", "control_test", "grid",
          "learner", "timeout")

        # submit to the cluster
        qsub_handle <- PRISM::qsub_lapply(
          X = qsub_x,
          qsub_environment = qsub_environment,
          qsub_packages = qsub_packages,
          qsub_config = qsub_config,
          FUN = qsub_fun
        )

        # save data and handle to RDS file
        out <- lst(
          method, obj_fun, design, tasks, task_group,
          task_fold, num_cores, metrics,
          control_train, control_test, grid, qsub_handle)
        saveRDS(out, qsubhandle_file)
      } else {
        # run locally
        out <- pbapply::pblapply(
          X = qsub_x,
          cl = num_cores,
          FUN = qsub_fun
        )
        out
      }
    }
  })
}

#' Downloading and processing the results of the benchmark jobs
#'
#' @param out_dir The folder in which to output intermediate and final results.
#'
#' @importFrom PRISM qsub_retrieve qacct qstat_j
#' @export
benchmark_suite_retrieve <- function(out_dir) {
  method_names <- list.dirs(out_dir, full.names = FALSE, recursive = FALSE) %>% discard(~ . == "")
  map_df(method_names, function(method_name) {
    method_folder <- paste0(out_dir, method_name)
    output_file <- paste0(method_folder, "/output.rds")
    qsubhandle_file <- paste0(method_folder, "/qsubhandle.rds")

    if (file.exists(output_file)) {
      cat(method_name, ": Reading previous output\n", sep = "")
      readRDS(output_file)
    } else if (!file.exists(output_file) && file.exists(qsubhandle_file)) {
      cat(method_name, ": Attempting to retrieve output from cluster: ", sep = "")
      data <- readRDS(qsubhandle_file)
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
        saveRDS(out_rds, output_file)
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
                 time_total = test_out$opt.path$env$exec.time,
                 param_i = seq_along(time_total)),
      test_out$opt.path$env$path
    ) %>% filter(param_i <= nrow(train_out$opt.path$env$path)) %>%
      mutate(iteration_i = train_out$opt.path$env$dob)
  )
  path_summary$params <- map(seq_len(nrow(path_summary)), ~extract_row_to_list(path_summary[,x_names], .))
  path_summary <- path_summary %>% select(-one_of(x_names))

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

