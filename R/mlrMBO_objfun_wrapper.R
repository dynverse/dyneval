#' Used for wrapping an evaluation function around a TI method
#'
#' @inheritParams execute_evaluation
#' @param noisy whether or not the metric is noisy or not
#'
#' @importFrom smoof makeSingleObjectiveFunction makeMultiObjectiveFunction
#' @export
make_obj_fun <- function(method, metrics, extra_metrics, noisy = FALSE) {
  # Use different makefunction if there are multiple metrics versus one
  if (length(metrics) > 1) {
    make_fun <- function(...) makeMultiObjectiveFunction(..., n.objectives = length(metrics))
  } else {
    make_fun <- makeSingleObjectiveFunction
  }

  # Wrap the method function in an evaluation function
  make_fun(
    name = "TItrain",
    vectorized = FALSE,
    minimize = rep(FALSE, length(metrics)),
    noisy = noisy,
    has.simple.signature = FALSE,
    par.set = method$par_set,
    fn = function(x, tasks, output_model, mc_cores = 1) {
      execute_evaluation(
        tasks = tasks,
        method = method,
        parameters = x,
        metrics = metrics,
        extra_metrics = extra_metrics,
        output_model = output_model,
        mc_cores = mc_cores
      )
    }
  )
}
