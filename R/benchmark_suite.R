#' @importFrom testthat expect_equal
#' @importFrom PRISM qsub_lapply override_qsub_config
#' @importFrom ParamHelpers generateDesignOfDefaults generateDesign
benchmark_suite_submit <- function(
  tasks, task_group, task_fold,
  methods = get_descriptions(as_tibble = F),
  num_cores = 4,
  metrics = c("auc_R_nx", "robbie_network_score"),
  num_iterations = 20,
  num_init_params = 100,
  output_root_folder = "results/"
  ) {
  testthat::expect_equal(nrow(tasks), length(task_group))
  testthat::expect_equal(nrow(tasks), length(task_fold))

  impute_fun <- impute_y_fun(length(metrics))

  ## MBO settings
  control_train <- makeMBOControl(
    n.objectives = length(metrics),
    propose.points = num_cores,
    impute.y.fun = impute_fun) %>%
    setMBOControlTermination(iters = num_iterations) %>%
    setMBOControlInfill(makeMBOInfillCritDIB())
  control_test <- control_train
  control_test$iters <- 1
  control_test$propose.points <- 1

  grid <- expand.grid(
    fold_i = unique(tasks$subtask_ix),
    group_sel = unique(select_tasks$group),
    repeat_i = seq_len(1),
    stringsAsFactors = F
  )

  ## Run MBO
  for (method in methods) {
    method_folder <- paste0(output_root_folder, method$short_name)
    output_file <- paste0(method_folder, "/output.rds")
    qsubhandle_file <- paste0(method_folder, "/qsubhandle.rds")

    dir.create(method_folder, recursive = T, showWarnings = F)

    if (!file.exists(output_file) && !file.exists(qsubhandle_file)) {
      cat("Submitting ", method$name, "\n", sep="")

      obj_fun <- make_obj_fun(method, metrics = metrics)

      design <- bind_rows(
        ParamHelpers::generateDesignOfDefaults(method$par_set),
        ParamHelpers::generateDesign(n = num_init_params, par.set = method$par_set)
      )

      qsub_handle <- PRISM::qsub_lapply(
        X = seq_len(nrow(grid)),
        qsub_environment = c(
          "method", "obj_fun", "design",
          "tasks", "task_group", "task_fold",
          "num_cores", "metrics", "impute_fun",
          "control_train", "control_test", "grid"),
        qsub_packages = c("dplyr", "purr", "dyneval", "mlrMBO", "parallelMap"),
        qsub_config = PRISM::override_qsub_config(
          wait = F,
          num_cores = num_cores,
          memory = "20G",
          name = paste0("D_", method$short_name),
          remove_tmp_folder = F,
          stop_on_error = F,
          verbose = F,
          max_wall_time = "99:00:00",
          execute_before = " export R_MAX_NUM_DLLS=200"
        ),
        FUN = function(grid_i) {
          fold_i <- grid[grid_i,]$fold_i
          group_sel <- grid[grid_i,]$group_sel
          repeat_i <- grid[grid_i,]$repeat_i

          ## start parameter optimisation
          parallelStartMulticore(cpus = num_cores, show.info = TRUE)
          tune_train <- mbo(
            obj_fun,
            design = design,
            control = control_train,
            show.info = T,
            more.args = list(tasks = select_tasks %>% filter(group == group_sel, subtask_ix != fold_i))
          )
          tune_test <- mbo(
            obj_fun,
            design = tune_train$opt.path$env$path %>% select(-starts_with("y_"), -one_of("y")),
            control = control_test,
            show.info = T,
            more.args = list(tasks = select_tasks %>% filter(group == group_sel, subtask_ix == fold_i))
          )
          parallelStop()

          list(design = design, tune_train = tune_train, tune_test = tune_test)
        })

      out <- lst(
        method, obj_fun, design, tasks, task_group,
        task_fold, num_cores, metrics, impute_fun,
        control_train, control_test, grid, qsub_handle)
      saveRDS(out, qsubhandle_file)
    }
  }
}

# benchmark_suite_retrieve <- function() {
#   ## Post process fun
#   post_fun <- function(rds_i, out_rds) {
#     grid_i <- rds_i
#     fold_i <- grid$fold_i[[rds_i]]
#     group_sel <- grid$group_sel[[rds_i]]
#     repeat_i <- grid$repeat_i[[rds_i]]
#     design <- out_rds$design
#     train_out <- out_rds$tune_train
#     test_out <- out_rds$tune_test
#
#     eval_summ_gath <- bind_rows(
#       data.frame(type = "train", train_out$opt.path$env$path) %>% mutate(
#         grid_i, repeat_i, fold_i, group_sel,
#         param_i = seq_len(n()),
#         time = train_out$opt.path$env$exec.time
#       ),
#       data.frame(type = "test", test_out$opt.path$env$path) %>% mutate(
#         grid_i, repeat_i, fold_i, group_sel,
#         param_i = seq_len(n()),
#         time = test_out$opt.path$env$exec.time
#       )
#     ) %>% filter(param_i <= nrow(train_out$opt.path$env$path)) %>%
#       dplyr::as_data_frame()
#
#     if (!all(eval_summ_gath$y_1 == -1)) {
#       eval_summ <- eval_summ_gath %>%
#         gather(eval_metric, score, y_1:y_3, time) %>%
#         mutate(comb = paste0(type, "_", eval_metric)) %>%
#         select(-type, -eval_metric) %>%
#         spread(comb, score) %>%
#         mutate(iteration_i = train_out$opt.path$env$dob[param_i]) %>%
#         arrange(param_i)
#
#       ## collect the scores per dataset individually
#       eval_ind <- bind_rows(lapply(seq_len(nrow(eval_summ)), function(param_i) {
#         iteration_i <- eval_summ$iteration_i[[param_i]]
#         bind_rows(
#           if (eval_summ$train_y_1[[param_i]] >= 0) { # did this execution finish correctly?
#             train_out$opt.path$env$extra[[param_i]]$.summary %>% mutate(grid_i, repeat_i, fold_i, group_sel, param_i, iteration_i, fold_type = "train")
#           } else {
#             NULL
#           },
#           if (eval_summ$test_y_1[[param_i]] >= 0) { # did this execution finish correctly?
#             test_out$opt.path$env$extra[[param_i]]$.summary %>% mutate(grid_i, repeat_i, fold_i, group_sel, param_i, iteration_i, fold_type = "test")
#           } else {
#             NULL
#           }
#         )
#       })) %>% left_join(tasks %>% dplyr::select(type, ti_type, id, experiment_id, platform_id, version, task_ix, subtask_ix), by = c("task_id" = "id")) %>%
#         as_data_frame
#
#       ## group them together per ti_type
#       eval_grp <- eval_ind %>% group_by(ti_type, grid_i, repeat_i, fold_i, group_sel, iteration_i, param_i, fold_type) %>% summarise_at(metrics, mean) %>% ungroup()
#     } else {
#       eval_summ <- NULL
#       eval_summ_gath <- NULL
#       eval_ind <- NULL
#       eval_grp <- NULL
#     }
#
#     dplyr::lst(eval_summ, eval_summ_gath, eval_ind, eval_grp)
#   }
#
#   ## Process results
#   for (method in methods) {
#     method_folder <- paste0(output_root_folder, method$short_name)
#     output_file <- paste0(method_folder, "/output.rds")
#     qsubhandle_file <- paste0(method_folder, "/qsubhandle.rds")
#
#     dir.create(method_folder, showWarnings = F)
#
#     if (!file.exists(output_file) && file.exists(qsubhandle_file)) {
#       qsub_handle <- readRDS(qsubhandle_file)
#
#       load(qsub_handle$src_rdata)
#
#       output <- qsub_retrieve(qsub_handle, post_fun = post_fun, wait = F)
#     }
#   }
# }
