
#' Compute metrics related to the prediction of the positions
#'
#' @param dataset A dataset containing a trajectory
#' @param prediction A predicted trajectory
#' @param metrics Which metrics to predict, can be rf_mse, rf_rsq, rf_nmse, lm_mse, lm_rsq and/or lm_nmse
#'
#' @keywords metric
#'
#' @importFrom reshape2 acast
#' @importFrom ranger ranger
#' @importFrom stats lm sd
calculate_position_predict <- function(dataset, prediction, metrics = c("rf_mse", "rf_rsq", "lm_mse", "lm_rsq")) {
  cell_ids <- dataset$cell_ids

  output <- list(
    summary = list()
  )

  gold_milenet_m <- dataset$milestone_percentages %>%
    reshape2::acast(cell_id ~ milestone_id, value.var = "percentage", fill = 0) %>%
    expand_matrix(rownames = cell_ids)

  # calculate the baseline mean squared error, by using the mean of each milestone as the prediction
  baseline_mse <- (t(gold_milenet_m) - apply(gold_milenet_m, 2, mean)) %>% apply(1, function(x) mean(x^2)) %>% mean()

  if (!is.null(prediction) && length(unique(prediction$milestone_percentages$cell_id)) >= 3) {
    pred_milenet_m <- prediction$milestone_percentages %>%
      reshape2::acast(cell_id ~ milestone_id, value.var = "percentage", fill = 0) %>%
      expand_matrix(rownames = cell_ids)
    pred_milenet_m <- pred_milenet_m[, apply(pred_milenet_m, 2, stats::sd) > 0]

    # random forest
    if (any(c("rf_mse", "rf_rsq", "rf_nmse") %in% metrics)) {
      rfs <- map(seq_len(ncol(gold_milenet_m)), function(i) {
        data <- gold_milenet_m[,i, drop = F] %>%
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

      output$rf_mses <- map_dbl(rfs, ~ mean(.$prediction.error)) %>% setNames(colnames(gold_milenet_m))
      output$summary$rf_mse <- mean(output$rf_mses)

      output$rf_rsqs <- map_dbl(rfs, ~ mean(.$r.squared)) %>% setNames(colnames(gold_milenet_m))
      output$rf_rsqs[is.na(output$rf_rsqs)] <- 1 # if no cells are nearby this milestones, the rsq will obviously be perfect
      output$summary$rf_rsq <- mean(output$rf_rsqs) %>% max(0)

      output$summary$rf_nmse <- (1 - output$summary$rf_mse / baseline_mse) %>% max(0)
    }

    # linear model
    if (any(c("lm_mse", "lm_rsq", "lm_nmse") %in% metrics)) {
      lms <- map(seq_len(ncol(gold_milenet_m)), function(i) {
        data <- gold_milenet_m[,i, drop = F] %>%
          magrittr::set_colnames("PREDICT") %>%
          cbind(pred_milenet_m) %>%
          as.data.frame()

        model <- stats::lm(PREDICT~., data=data)

        list(
          mse = mean(model$residuals^2),
          rsq = suppressWarnings(summary(model))$r.squared
        )
      })

      output$summary$lm_mse <- map_dbl(lms, "mse") %>% mean()

      output$lm_rsqs <- map_dbl(lms, "rsq")
      output$lm_rsqs[is.na(output$lm_rsqs)] <- 1 # if no cells are nearby this milestones, the rsq will obviously be perfect
      output$summary$lm_rsq <- mean(output$lm_rsqs) %>% max(0)

      output$summary$lm_nmse <- (1 - output$summary$lm_mse / baseline_mse) %>% max(0)
    }
  } else {
    output$summary <- list(
      rf_mse = baseline_mse,
      rf_nmse = 0,
      rf_rsq = 0,
      lm_mse = baseline_mse,
      lm_rsq = 0,
      lm_nmse = 0
    )
  }

  output
}
