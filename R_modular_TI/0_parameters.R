generate_parameters <- function(parameter_sets) {
  dplyr::bind_rows(lapply(parameter_sets, function(param_set) {
    do.call("expand.grid", c(param_set, list(stringsAsFactors = F)))
  }))
}

#' Return a set of parameter values from the expanded parameter grid
#'
#' @param method the wrapped method to extract the expanded parameters from
#' @param parameter_row the index of the expanded parameters to use
#'
#' @export
get_parameter_row <- function(method, parameter_row) {

  na.omit(as.list(method$expanded_parameters[parameter_row,,drop=F]))
}
