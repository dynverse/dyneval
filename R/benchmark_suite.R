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
#'
#' @importFrom testthat expect_equal
#' @importFrom PRISM qsub_lapply override_qsub_config
#' @importFrom mlrMBO makeMBOControl setMBOControlTermination setMBOControlInfill makeMBOInfillCritDIB
#' @importFrom mlr makeLearner
#' @importFrom ParamHelpers generateDesignOfDefaults generateDesign
#' @importFrom parallelMap parallelStartMulticore parallelStop
#' @importFrom randomForest randomForest
#'
#' @export
benchmark_suite_submit <- function(
  tasks,
  task_group,
  task_fold,
  out_dir,
  timeout = 60 * nrow(tasks),
  methods = get_descriptions(as_tibble = TRUE),
  metrics = c("auc_R_nx", "robbie_network_score"),
  num_cores = 4,
  memory = "20G",
  max_wall_time = "72:00:00",
  num_iterations = 20,
  num_init_params = 100,
  num_repeats = 1
) {
  testthat::expect_is(tasks, "tbl")
  testthat::expect_equal(nrow(tasks), length(task_group))
  testthat::expect_equal(nrow(tasks), length(task_fold))
  testthat::expect_is(methods, "tbl")

  ## MBO settings
  control_train <- mlrMBO::makeMBOControl(
    n.objectives = length(metrics),
    propose.points = num_cores,
    impute.y.fun = impute_y_fun(length(metrics))) %>%
    mlrMBO::setMBOControlTermination(iters = num_iterations) %>%
    mlrMBO::setMBOControlInfill(mlrMBO::makeMBOInfillCritDIB())
  control_test <- control_train
  control_test$iters <- 1
  control_test$propose.points <- 1
  learner <- mlr::makeLearner(
    "regr.randomForest",
    se.method = "jackknife",
    predict.type = "se",
    keep.inbag = TRUE)

  ## Grid settings
  grid <- expand.grid(
    fold_i = sort(unique(task_fold)),
    group_sel = sort(unique(task_group)),
    repeat_i = seq_len(num_repeats),
    stringsAsFactors = F
  )

  ## Run MBO
  for (methodi in seq_len(nrow(methods))) {
    method <- dynutils::extract_row_to_list(methods, methodi)

    # determine where to store certain outputs
    method_folder <- paste0(out_dir, method$short_name)
    output_file <- paste0(method_folder, "/output.rds")
    qsubhandle_file <- paste0(method_folder, "/qsubhandle.rds")

    dir.create(method_folder, recursive = TRUE, showWarnings = FALSE)

    if (!file.exists(output_file) && !file.exists(qsubhandle_file)) {
      cat("Submitting ", method$name, "\n", sep="")

      # create an objective function
      obj_fun <- make_obj_fun(method = method, metrics = metrics, timeout = timeout)

      # generate initial parameters
      design <- bind_rows(
        ParamHelpers::generateDesignOfDefaults(method$par_set),
        ParamHelpers::generateDesign(n = num_init_params, par.set = method$par_set)
      )

      # cluster params
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
      qsub_packages <- c("dplyr", "purrr", "dyneval", "mlrMBO", "parallelMap")
      qsub_environment <-  c(
        "method", "obj_fun", "design",
        "tasks", "task_group", "task_fold",
        "num_cores", "metrics",
        "control_train", "control_test", "grid")

      # submit to the cluster
      qsub_handle <- PRISM::qsub_lapply(
        X = seq_len(nrow(grid)),
        qsub_environment = qsub_environment,
        qsub_packages = qsub_packages,
        qsub_config = qsub_config,
        FUN = function(grid_i) {
          fold_i <- grid[grid_i,]$fold_i
          group_sel <- grid[grid_i,]$group_sel
          repeat_i <- grid[grid_i,]$repeat_i

          ## start parameter optimisation
          parallelStartMulticore(cpus = num_cores, show.info = TRUE)
          tune_train <- mbo(
            obj_fun,
            learner = learner,
            design = design,
            control = control_train,
            show.info = TRUE,
            more.args = list(tasks = tasks[task_group == group_sel & task_fold != fold_i,])
          )
          tune_test <- mbo(
            obj_fun,
            learner = learner,
            design = tune_train$opt.path$env$path %>% select(-starts_with("y_"), -one_of("y")),
            control = control_test,
            show.info = TRUE,
            more.args = list(tasks = tasks[task_group == group_sel & task_fold == fold_i,])
          )
          parallelStop()

          list(design = design, tune_train = tune_train, tune_test = tune_test)
        })

      out <- lst(
        method, obj_fun, design, tasks, task_group,
        task_fold, num_cores, metrics,
        control_train, control_test, grid, qsub_handle)
      saveRDS(out, qsubhandle_file)
    }
  }
}

#' Downloading and processing the results of the benchmark jobs
#'
#' @param out_dir The folder in which to output intermediate and final results.
#'
#' @importFrom PRISM qsub_retrieve qacct
#' @export
benchmark_suite_retrieve <- function(out_dir) {
  method_names <- list.files(out_dir)
  for (method_name in method_names) {
    method_folder <- paste0(out_dir, method_name)
    output_file <- paste0(method_folder, "/output.rds")
    qsubhandle_file <- paste0(method_folder, "/qsubhandle.rds")

    if (!file.exists(output_file) && file.exists(qsubhandle_file)) {
      data <- readRDS(qsubhandle_file)
      qsub_handle <- data$qsub_handle

      output <- qsub_retrieve(
        qsub_handle,
        post_fun = function(rds_i, out_rds) {
          benchmark_suite_retrieve_helper(rds_i, out_rds, data)
        },
        wait = FALSE
      )

      if (!is.null(output)) {
        cat("Saving output of ", method_name, "\n", sep = "")

        which_errored <- sapply(output, is.na)
        outputs <- lapply(output, function(x) if(is.na(x)) NULL else x)
        errors <- lapply(output, function(x) if(!is.na(x)) NULL else attr(x, "qsub_error"))
        qacct <- qacct(qsub_handle)

        rds_lst <- lst(
          outputs,
          which_errored,
          errors,
          qacct,
          qsub_handle
        )
        saveRDS(output, output_file)
      }
    }
  }
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

  eval_summ_gath <- bind_rows(
    data.frame(type = "train", train_out$opt.path$env$path, stringsAsFactors = F) %>%
      mutate(
      grid_i, repeat_i, fold_i, group_sel,
      param_i = seq_len(n()),
      time = train_out$opt.path$env$exec.time
    ),
    data.frame(type = "test", test_out$opt.path$env$path, stringsAsFactors = F) %>%
      mutate(
      grid_i, repeat_i, fold_i, group_sel,
      param_i = seq_len(n()),
      time = test_out$opt.path$env$exec.time
    )
  ) %>% filter(param_i <= nrow(train_out$opt.path$env$path)) %>%
    as_data_frame()

  if (!all(eval_summ_gath$y_1 == -1)) {
    eval_summ <- eval_summ_gath %>%
      gather(eval_metric, score, starts_with("y_"), starts_with("time")) %>%
      mutate(comb = paste0(type, "_", eval_metric)) %>%
      select(-type, -eval_metric) %>%
      spread(comb, score) %>%
      mutate(iteration_i = train_out$opt.path$env$dob[param_i]) %>%
      arrange(param_i)

    ## collect the scores per dataset individually
    eval_ind <- bind_rows(lapply(seq_len(nrow(eval_summ)), function(param_i) {
      iteration_i <- eval_summ$iteration_i[[param_i]]
      bind_rows(
        if (eval_summ$train_y_1[[param_i]] >= 0) { # did this execution finish correctly?
          train_out$opt.path$env$extra[[param_i]]$.summary %>%
            mutate(grid_i, repeat_i, fold_i, group_sel, param_i, iteration_i, fold_type = "train")
        } else {
          NULL
        },
        if (eval_summ$test_y_1[[param_i]] >= 0) { # did this execution finish correctly?
          test_out$opt.path$env$extra[[param_i]]$.summary %>%
            mutate(grid_i, repeat_i, fold_i, group_sel, param_i, iteration_i, fold_type = "test")
        } else {
          NULL
        }
      )
    })) %>% as_data_frame

    ## group them together per task_group
    eval_grp <- eval_ind %>%
      left_join(tasks %>% mutate(task_fold, task_group) %>% select(task_id = id, ti_type, task_fold, task_group), by = "task_id") %>%
      group_by(
        ti_type, task_group, grid_i, repeat_i, fold_i,
        group_sel, iteration_i, param_i, fold_type,
        method_name, method_short_name
      ) %>%
      select(-task_id, -task_fold) %>%
      summarise_all(mean) %>% ungroup()

    lst(eval_summ, eval_summ_gath, eval_ind, eval_grp)
  } else {
    x <- NA
    attr(x, "qsub_error") <- "All parameters produced errors, resulting in a default error score"
    x
  }


}
