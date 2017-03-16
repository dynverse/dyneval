imp_similarity_to_distance <- list(
  wrap_method(
    method_name = "smacof_sim2diss",
    method_type = mt_similarity_to_distance,
    method_function = function(similarity, method) {
      list(
        distance = smacof::sim2diss(similarity, method = method)
      )
    },
    parameter_sets = list(list(method = c("corr", "neglog", "counts"))),
    required_namespaces = c("smacof")
  )
)
imp_symmetric_similarity_to_distance <- list(
  wrap_method(
    method_name = "smacof_sim2diss",
    method_type = mt_symmetric_similarity_to_distance,
    method_function = function(similarity, method) {
      list(
        distance = smacof::sim2diss(similarity, method = method)
      )
    },
    parameter_sets = list(list(method = c("corr", "neglog", "counts"))),
    required_namespaces = c("smacof")
  )
)
