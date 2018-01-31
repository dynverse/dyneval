
#' Compute the random forest OOB-MSE metric
#'
#' @param task A task
#' @param prediction A predicted model
#'
#' @importFrom randomForest randomForest
#' @importFrom reshape2 acast
#' @importFrom tibble lst
compute_rfmse <- function(task, prediction) {
  cell_ids <- task$cell_ids

  if (!is.null(prediction)) {
    gold_milenet_m <- task$milestone_percentages %>%
      reshape2::acast(cell_id ~ milestone_id, value.var = "percentage", fill = 0) %>%
      expand_mat(cell_ids)
    pred_milenet_m <- prediction$milestone_percentages %>%
      reshape2::acast(cell_id ~ milestone_id, value.var = "percentage", fill = 0) %>%
      expand_mat(cell_ids)

    rfs <- lapply(seq_len(ncol(gold_milenet_m)), function(i) {
      randomForest::randomForest(pred_milenet_m, gold_milenet_m[,i])
    })

    mses <- map_dbl(rfs, ~ mean(.$mse)) %>% setNames(colnames(gold_milenet_m))
    mmse <- mean(mses)

    rsqs <- map_dbl(rfs, ~ mean(.$rsq)) %>% setNames(colnames(gold_milenet_m))
    mrsq <- mean(rsqs)

    summary <- lst(rf_mse, rf_rsq)

    lst(
      mses,
      rsqs,
      summary
    )
  } else {
    lst(
      summary = lst(
        rf_mse = Inf,
        rf_rsq = 0
      )
    )
  }
}

expand_mat <- function(mat, rownames) {
  newmat <- matrix(0, nrow = length(rownames), ncol = ncol(mat), dimnames = list(rownames, colnames(mat)))
  newmat[rownames(mat),] <- mat
  newmat
}
