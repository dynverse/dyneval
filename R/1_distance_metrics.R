imp_symmetric_distance_metric <- list(
  wrap_method(
    method_name = "symmetric_euclidean",
    method_type = mt_symmetric_distance_metric,
    method_function = function(x) {
      dm <- as.matrix(stats::dist(x))
      list(
        distance = wrap_data_object(dt_symmetric_distance_matrix, dm)
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
      dm <- SCORPIUS:::euclidean_distance_rcpp(x, y)
      list(
        distance = wrap_data_object(dt_distance_matrix, dm)
      )
    },
    parameter_sets = list(),
    required_namespaces = c("SCORPIUS")
  )
)

# simlr = function(x, ndim=3, nclusters=4) {
#   result = SIMLR::SIMLR(t(x), nclusters)
#   S = result$S
#   space = tsne::tsne(as.dist(max(S)-S), k = ndim)
#   process_dimred(space)
# }
