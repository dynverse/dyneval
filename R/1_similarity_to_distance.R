imp_similarity_to_distance <- list(
  wrap_method(
    method_name = "smacof_sim2diss",
    method_type = mt_similarity_to_distance,
    method_function = function(similarity, method) {
      dm <- smacof::sim2diss(similarity, method = method)
      list(
        distance = wrap_data_object(dt_distance_matrix, dm)
      )
    },
    parameter_sets = list(list(method = c("corr", "neglog", "counts"))),
    required_namespaces = c("smacof")
  )
)
