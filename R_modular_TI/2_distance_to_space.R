imp_distance_to_space <- list(
  wrap_method(
    method_name = "mds_stats",
    method_type = "distance_dimensionality_reduction",
    method_function = function(distance, num_dimensions) {
      sp <- stats::cmdscale(distance, k = num_dimensions)
      list(
        space = clean_dimred_output(sp, rownames(distance))
      )
    },
    parameter_sets = list(list(num_dimensions = seq(2, 10))),
    required_namespaces = c("stats")
  ),
  wrap_method(
    method_name = "mds_sammon",
    method_type = "distance_dimensionality_reduction",
    method_function = function(distance, num_dimensions) {
      sp <- MASS::sammon(distance, k = num_dimensions)$points
      list(
        space = clean_dimred_output(sp, rownames(distance))
      )
    },
    parameter_sets = list(list(num_dimensions = seq(2, 10))),
    required_namespaces = c("MASS")
  ),
  wrap_method(
    method_name = "mds_iso",
    method_type = "distance_dimensionality_reduction",
    method_function = function(distance, num_dimensions) {
      sp <- MASS::isoMDS(distance, k = num_dimensions)$points
      list(
        space = clean_dimred_output(sp, rownames(distance))
      )
    },
    parameter_sets = list(list(num_dimensions = seq(2, 10))),
    required_namespaces = c("MASS")
  ),
  wrap_method(
    method_name = "mds_smacof",
    method_type = "distance_dimensionality_reduction",
    method_function = function(distance, num_dimensions) {
      sp <- smacof::mds(distance, ndim = num_dimensions)$conf
      list(
        space = clean_dimred_output(sp, rownames(distance))
      )
    },
    parameter_sets = list(list(num_dimensions = seq(2, 10))),
    required_namespaces = c("smacof")
  ),
  wrap_method(
    method_name = "tsne_tsne",
    method_type = "distance_dimensionality_reduction",
    method_function = function(distance, num_dimensions) {
      sp <- tsne::tsne(as.dist(distance), k = num_dimensions)
      list(
        space = clean_dimred_output(sp, rownames(distance))
      )
    },
    parameter_sets = list(list(num_dimensions = seq(2, 10))),
    required_namespaces = c("tsne")
  ),
  wrap_method(
    method_name = "Rtsne_Rtsne",
    method_type = "distance_dimensionality_reduction",
    method_function = function(distance, num_dimensions) {
      sp <- Rtsne::Rtsne(as.dist(distance), dims = num_dimensions, is_distance = T)
      list(
        space = clean_dimred_output(sp, rownames(distance))
      )
    },
    parameter_sets = list(list(num_dimensions = seq(2, 10))),
    required_namespaces = c("Rtsne")
  ),
  wrap_method(
    method_name = "diffusionmap_diffusionMap",
    method_type = "distance_dimensionality_reduction",
    method_function = function(distance, num_dimensions, num_eigen) {
      sp <- diffusionMap::diffuse(as.dist(distance), neigen = num_eigen)$X[,seq_len(num_dimensions)]
      list(
        space = clean_dimred_output(sp, rownames(distance))
      )
    },
    parameter_sets = list(list(num_dimensions = seq(2, 10), num_eigen = NULL)),
    required_namespaces = c("diffusionMap")
  )
)
