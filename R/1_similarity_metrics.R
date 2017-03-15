imp_symmetric_similarity_metric <- list(
  wrap_method(
    method_name = "symmetric_correlation",
    method_type = mt_symmetric_similarity_metric,
    method_function = function(x, method) {
      sm <- stats::cor(t(x), method = method)
      list(
        similarity = wrap_data_object(dt_symmetric_similarity_matrix, sm)
      )
    },
    parameter_sets = list(list(method = c("pearson", "kendall", "spearman"))),
    required_namespaces = c("stats")
  )
)

imp_similarity_metric <- list(
  wrap_method(
    method_name = "correlation",
    method_type = mt_similarity_metric,
    method_function = function(x, y, method) {
      sm <- stats::cor(t(x), t(y), method = method)
      list(
        similarity = wrap_data_object(dt_similarity_matrix, sm)
      )
    },
    parameter_sets = list(list(method = c("pearson", "kendall", "spearman"))),
    required_namespaces = c("stats")
  )
)

