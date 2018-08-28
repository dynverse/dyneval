
#' Compare feature importances derived by both trajectories
#'
#' @param dataset A dataset
#' @param prediction A predicted model
#' @param num_trees The number of trees to use for the random forest
#' @param mtry Number of features to split in each node. Can be a function with as argument the dataset
#'
#' @importFrom dynfeature calculate_overall_feature_importance
#'
#' @export
calculate_featureimp_cor <- function(dataset, prediction, num_trees = 10000, mtry = function(x) ncol(x) * .01) {
  cell_ids <- dataset$cell_ids

  if (!is.null(prediction) && length(unique(prediction$milestone_percentages$cell_id)) >= 3) {
    method_params <- list(num.trees = num_trees, mtry = mtry)

    dataset_imp <- dynfeature::calculate_overall_feature_importance(dataset, expression_source = dataset$expression, method_params = method_params)
    pred_imp <- dynfeature::calculate_overall_feature_importance(prediction, expression_source = dataset$expression, method_params = method_params)

    join <- full_join(
      dataset_imp %>% rename(dataset_imp = importance),
      pred_imp %>% rename(pred_imp = importance),
      by = "feature_id"
    ) %>%
      mutate_at(c("dataset_imp", "pred_imp"), ~ ifelse(is.na(.), 0, .))

    featureimp_cor <- cor(join$dataset_imp, join$pred_imp) %>% max(0)

    cov_wt <- cov.wt(
      x = matrix(c(join$dataset_imp, join$pred_imp), ncol = 2),
      wt = join$dataset_imp,
      cor = TRUE
    )
    featureimp_wcor <- cov_wt$cor[1, 2] %>% max(0)

    lst(
      featureimp_cor,
      featureimp_wcor
    )
  } else {
    list(
      featureimp_cor = 0,
      featureimp_wcor = 0
    )
  }
}


#' Compare enrichment in finding back the most important genes
#'
#' @param dataset A dataset
#' @param prediction A predicted model
#' @param num_trees the number of trees to use during the calculation of the metric
#'
#' @importFrom dynfeature calculate_overall_feature_importance
#' @importFrom stats ks.test wilcox.test
calculate_featureimp_enrichment <- function(dataset, prediction, num_trees = 10000) {
  cell_ids <- dataset$cell_ids

  tryCatch({
    if (!is.null(prediction) && length(unique(prediction$milestone_percentages$cell_id)) >= 3) {
      method_params <- list(num.trees = num_trees)
      pred_imp <- dynfeature::calculate_overall_feature_importance(prediction, expression_source = dataset$expression, method_params = method_params)
      dataset_features <- dataset$prior_information$features_id

      sel <- pred_imp$importance[pred_imp$feature_id %in% dataset_features]
      notsel <- pred_imp$importance[!pred_imp$feature_id %in% dataset_features]

      if (length(notsel) > 2) {
        ks <- stats::ks.test(sel, notsel, alternative = "greater")
        wilcox <- stats::wilcox.test(sel, notsel, alternative = "greater")

        list(
          featureimp_ks = ks$p.value,
          featureimp_wilcox = 1 - wilcox$p.value
        )
      } else {
        list(featureimp_ks = 1, featureimp_wilcox = 1)
      }
    } else {
      list(featureimp_ks = 0, featureimp_wilcox = 0)
    }
  }, error = function(e) {
    warning("featureimp_enrichment errored! check reason!")
    list(featureimp_ks = 0, featureimp_wilcox = 0)
  })
}

#' @examples
#' dataset <- dyntoy::generate_dataset(num_cells = 300, num_features = 300)
#' prediction <- dynwrap::infer_trajectory(dataset, "slingshot", parameters = list())
#' num_trees <- 10000
#' mtry = function(x) ncol(x) * .01
#' calculate_featureimp_cor(dataset, prediction)
