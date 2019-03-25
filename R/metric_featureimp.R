
#' Compare feature importances derived by both trajectories
#'
#' @param dataset A dataset
#' @param prediction A predicted trajectory
#' @param expression_source The expression data matrix, with features as columns.
#'   * If a matrix is provided, it is used as is.
#'   * If a character is provided, `dataset[[expression_source]]` should contain the matrix.
#'   * If a function is provided, that function will be called in order to obtain the expression (useful for lazy loading).
#'
#' @keywords metric
#'
#' @inheritParams dynfeature::calculate_overall_feature_importance
#'
#' @importFrom dynfeature calculate_overall_feature_importance fi_ranger_rf_lite
#'
#' @export
calculate_featureimp_cor <- function(
  dataset,
  prediction,
  expression_source = dataset$expression_source,
  fi_method = dynfeature::fi_ranger_rf_lite()
) {
  if (!is.null(prediction) && length(unique(prediction$milestone_percentages$cell_id)) >= 3) {
    dataset_imp <-
      dynfeature::calculate_overall_feature_importance(
        trajectory = dataset,
        expression_source = expression_source,
        fi_method = fi_method
      )
    pred_imp <-
      dynfeature::calculate_overall_feature_importance(
        trajectory = prediction,
        expression_source = expression_source,
        fi_method = fi_method
      )

    .calculate_featureimp_cor(dataset_imp, pred_imp)
  } else {
    list(
      featureimp_cor = 0,
      featureimp_wcor = 0
    )
  }
}

.calculate_featureimp_cor <- function(dataset_imp, pred_imp) {
  join <- full_join(
    dataset_imp %>% rename(dataset_imp = importance),
    pred_imp %>% rename(pred_imp = importance),
    by = "feature_id"
  ) %>%
    mutate_at(c("dataset_imp", "pred_imp"), ~ ifelse(is.na(.), 0, .))

  # if either of the feauture importances are all the same, return 0 as correlation
  if (sd(join$dataset_imp) == 0 || sd(join$pred_imp) == 0) {
    list(
      featureimp_cor = 0,
      featureimp_wcor = 0
    )
  } else {
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
  }
}

#' Compare enrichment in finding back the most important genes
#'
#' @param dataset A dataset
#' @param prediction A predicted trajectory
#' @param expression_source The expression data matrix, with features as columns.
#'   * If a matrix is provided, it is used as is.
#'   * If a character is provided, `dataset[[expression_source]]` should contain the matrix.
#'   * If a function is provided, that function will be called in order to obtain the expression (useful for lazy loading).
#' @inheritParams dynfeature::calculate_overall_feature_importance
#'
#' @keywords metric
#'
#' @importFrom dynfeature calculate_overall_feature_importance fi_ranger_rf_lite
#' @importFrom stats ks.test wilcox.test
calculate_featureimp_enrichment <- function(
  dataset,
  prediction,
  expression_source = dataset$expression,
  fi_method = dynfeature::fi_ranger_rf_lite()
) {
  tryCatch({
    if (!is.null(prediction) && length(unique(prediction$milestone_percentages$cell_id)) >= 3) {
      pred_imp <-
        dynfeature::calculate_overall_feature_importance(
          prediction,
          expression_source = expression_source,
          fi_method = fi_method
        )

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
