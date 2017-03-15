clean_dimred_output <- function(space, rownames) {
  space <- as.matrix(space)
  dimnames(space) <- list(rownames, paste0("Comp", seq_len(ncol(space))))
  space
}

imp_distance_to_space <- list(
  wrap_method(
    method_name = "mds_stats",
    method_type = mt_distance_dimensionality_reduction,
    method_function = function(distance, num_dimensions) {
      sp <- stats::cmdscale(distance, k = num_dimensions)
      spc <- clean_dimred_output(sp, rownames(distance))
      list(
        space = wrap_data_object(dt_reduced_space, spc)
      )
    },
    parameter_sets = list(list(num_dimensions = seq(2, 10))),
    required_namespaces = c("stats")
  ),
  wrap_method(
    method_name = "mds_sammon",
    method_type = mt_distance_dimensionality_reduction,
    method_function = function(distance, num_dimensions) {
      sp <- MASS::sammon(distance, k = num_dimensions)$points
      spc <- clean_dimred_output(sp, rownames(distance))
      list(
        space = wrap_data_object(dt_reduced_space, spc)
      )
    },
    parameter_sets = list(list(num_dimensions = seq(2, 10))),
    required_namespaces = c("MASS")
  ),
  wrap_method(
    method_name = "mds_iso",
    method_type = mt_distance_dimensionality_reduction,
    method_function = function(distance, num_dimensions) {
      sp <- MASS::isoMDS(distance, k = num_dimensions)$points
      spc <- clean_dimred_output(sp, rownames(distance))
      list(
        space = wrap_data_object(dt_reduced_space, spc)
      )
    },
    parameter_sets = list(list(num_dimensions = seq(2, 10))),
    required_namespaces = c("MASS")
  ),
  wrap_method(
    method_name = "mds_smacof",
    method_type = mt_distance_dimensionality_reduction,
    method_function = function(distance, num_dimensions) {
      sp <- smacof::mds(distance, k = num_dimensions)$points
      spc <- clean_dimred_output(sp, rownames(distance))
      list(
        space = wrap_data_object(dt_reduced_space, spc)
      )
    },
    parameter_sets = list(list(num_dimensions = seq(2, 10))),
    required_namespaces = c("smacof")
  ),
  wrap_method(
    method_name = "tsne_tsne",
    method_type = mt_distance_dimensionality_reduction,
    method_function = function(distance, num_dimensions) {
      sp <- tsne::tsne(as.dist(distance), k = num_dimensions)
      spc <- clean_dimred_output(sp, rownames(distance))
      list(
        space = wrap_data_object(dt_reduced_space, spc)
      )
    },
    parameter_sets = list(list(num_dimensions = seq(2, 10))),
    required_namespaces = c("tsne")
  ),
  wrap_method(
    method_name = "diffusionmap_diffusionMap",
    method_type = mt_distance_dimensionality_reduction,
    method_function = function(distance, num_dimensions, num_eigen) {
      sp <- diffusionMap::diffuse(as.dist(distance), neigen = num_eigen)$X[,seq_len(num_dimensions)]
      spc <- clean_dimred_output(sp, rownames(distance))
      list(
        space = wrap_data_object(dt_reduced_space, spc)
      )
    },
    parameter_sets = list(list(num_dimensions = seq(2, 10), num_eigen = NULL)),
    required_namespaces = c("diffusionMap")
  )
)
