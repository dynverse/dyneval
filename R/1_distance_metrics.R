imp_symmetric_distance_metric <- list(
  wrap_method(
    method_name = "symmetric_euclidean",
    method_type = mt_symmetric_distance_metric,
    method_function = function(x) {
      list(
        distance = as.matrix(stats::dist(x))
      )
    },
    parameter_sets = list(),
    required_namespaces = c("stats")
  )
)

imp_distance_metric <- list(
  wrap_method(
    method_name = "euclidean",
    method_type = mt_distance_metric,
    method_function = function(x, y) {
      list(
        distance = SCORPIUS:::euclidean_distance_rcpp(x, y)
      )
    },
    parameter_sets = list(),
    required_namespaces = c("SCORPIUS")
  )
)
