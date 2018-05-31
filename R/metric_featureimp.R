
#' Compare feature importances derived by both trajectories
#'
#' @param task A task
#' @param prediction A predicted model
#' @param num_trees the number of trees to use during the calculation of the metric
#'
#' @importFrom dynfeature calculate_overall_feature_importance
compute_featureimp <- function(task, prediction, num_trees = 50000) {
  cell_ids <- task$cell_ids

  if (!is.null(prediction)) {

    method_params <- list(num.trees = num_trees)
    task_imp <- dynfeature::calculate_overall_feature_importance(task, expression_source = task$expression, method_params = method_params)
    pred_imp <- dynfeature::calculate_overall_feature_importance(prediction, expression_source = task$expression, method_params = method_params)

    imp_joined <- full_join(
      task_imp %>% rename(task_imp = importance),
      pred_imp %>% rename(pred_imp = importance),
      by = "feature_id"
    ) %>%
      mutate_at(c("task_imp", "pred_imp"), ~ ifelse(is.na(.), 0, .))

    # plot:
    # ggplot(imp_joined) + geom_point(aes(task_imp, pred_imp))
    featureimp_cor <- cor(imp_joined$task_imp, imp_joined$pred_imp)

    lst(
      featureimp_cor
    )
  } else {
    list(
      featureimp_cor = 0
    )
  }
}
