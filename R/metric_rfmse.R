
#' Compute the random forest OOB-MSE metric
#'
#' @param task A task
#' @param prediction A predicted model
#'
#' @importFrom reshape2 acast
#' @importFrom ranger ranger
compute_rfmse <- function(task, prediction) {
  cell_ids <- task$cell_ids

  if (!is.null(prediction)) {

    gold_milenet_m <- task$milestone_percentages %>%
      reshape2::acast(cell_id ~ milestone_id, value.var = "percentage", fill = 0) %>%
      expand_matrix(rownames = cell_ids)
    pred_milenet_m <- prediction$milestone_percentages %>%
      reshape2::acast(cell_id ~ milestone_id, value.var = "percentage", fill = 0) %>%
      expand_matrix(rownames = cell_ids)

    rfs <- lapply(seq_len(ncol(gold_milenet_m)), function(i) {
      data <- gold_milenet_m[,i, drop=F] %>%
        magrittr::set_colnames("PREDICT") %>%
        cbind(pred_milenet_m) %>%
        as.data.frame()

      ranger::ranger(
        dependent.variable.name = "PREDICT",
        data = data,
        num.trees = 5000,
        num.threads = 1
      )
    })

    mses <- map_dbl(rfs, ~ mean(.$prediction.error)) %>% setNames(colnames(gold_milenet_m))
    rf_mse <- mean(mses)

    rsqs <- map_dbl(rfs, ~ mean(.$r.squared)) %>% setNames(colnames(gold_milenet_m))
    rf_rsq <- mean(rsqs)

    summary <- lst(rf_mse, rf_rsq)

    lst(
      mses,
      rsqs,
      summary
    )
  } else {
    lst(
      summary = lst(
        rf_mse = 1,
        rf_rsq = 0
      )
    )
  }
}
