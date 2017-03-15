imp_matrix_to_space <- list(
  wrap_method(
    method_name = "pca_stats",
    method_type = mt_matrix_dimensionality_reduction,
    method_function = function(x, num_dimensions) {
      sp <- stats::prcomp(t(x))$rotation[,seq_len(num_dimensions)]
      spc <- clean_dimred_output(sp, rownames(distance))
      list(
        space = wrap_data_object(dt_reduced_space, spc)
      )
    },
    parameter_sets = list(list(num_dimensions = seq(2, 10))),
    required_namespaces = c("stats")
  ),
  wrap_method(
    method_name = "ica_fastICA",
    method_type = mt_matrix_dimensionality_reduction,
    method_function = function(x, num_dimensions) {
      x <- t(scale(t(x)))
      sp <- fastICA::fastICA(x, n.comp = num_dimensions)$S
      spc <- clean_dimred_output(sp, rownames(distance))
      list(
        space = wrap_data_object(dt_reduced_space, spc)
      )
    },
    parameter_sets = list(list(num_dimensions = seq(2, 10))),
    required_namespaces = c("fastICA")
  ),
  wrap_method(
    method_name = "lle_lle",
    method_type = mt_matrix_dimensionality_reduction,
    method_function = function(x, num_dimensions) {
      x <- t(scale(t(x)))
      num_neighbours <- lle::calc_k(t(scale(t(x))), num_dimensions)
      num_neighbours <- num_neighbours$k[which.min(num_neighbours$rho)]
      sp <- lle::lle(x, m = num_dimensions, k = num_neighbours)$Y
      spc <- clean_dimred_output(sp, rownames(distance))
      list(
        space = wrap_data_object(dt_reduced_space, spc)
      )
    },
    parameter_sets = list(list(num_dimensions = seq(2, 10))),
    required_namespaces = c("lle")
  )
)
