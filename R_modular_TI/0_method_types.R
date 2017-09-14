create_method_type <- function(name, input_types, output_types) {
  for (input_type in input_types) {
    check_data_type(input_type)
  }
  for (output_type in output_types) {
    check_data_type(output_type)
  }
  l <- list(
    name = name,
    input_types = input_types,
    output_types = output_types
  )
  class(l) <- "dyneval::method_type"
  l
}

wrap_method <- function(method_name, method_type, method_function, parameter_sets, required_namespaces) {
  check_method_type(method_type)
  method_type_object <- mt_objects[[method_type]]
  input_types <- method_type_object$input_types
  output_types <- method_type_object$output_types
  form_args <- formalArgs(method_function)
  if (!all(names(input_types) == form_args[seq_along(input_types)])) {
    stop("Not all input_types were found")
  }
  wrapped_function <- function(data_objects, extra_parameters) {
    for (namespace in required_namespaces) {
      requireNamespace(namespace)
    }
    for (i in seq_along(input_types)) {
      if (!instanceof(data_objects[[i]], input_types[[i]])) {
        message <- paste0(
          "Mismatch in input data types.\n",
          "Parameter: ", names(data_objects)[[i]], "\n",
          "Expected: ", paste(input_types[[i]], collapse=";"), "\n",
          "Observed: ", paste(data_objects[[i]], collapse=";"))
        stop(message)
      }
    }
    # call function with auto unwrap
    method_output <- do.call(method_function, c(lapply(data_objects, unwrap_data_object), extra_parameters))
    if (!all(names(output_types) == names(method_output))) {
      message <- paste0(
        "Not all output_types were found.\n",
        "Expected: ", paste(names(output_types), collapse = ", "), "\n",
        "Observed: ", paste(names(method_output), collapse = ", ")
      )
      stop(message)
    }
    # return with auto wrap
    # wrapped_output <- list()
    # for (i in seq_along(method_output)) {
    #   wrapped_output[[i]] <- wrap_data_object(output_types[[i]], method_output[[i]])
    # }
    # wrapped_output
    mapply(output_types, method_output, FUN = wrap_data_object, SIMPLIFY = FALSE)
  }
  expanded_parameters <- generate_parameters(parameter_sets)
  wrapped_method <- list(
    method_type = method_type,
    method_name = method_name,
    method_function = wrapped_function,
    parameter_sets = parameter_sets,
    expanded_parameters = expanded_parameters,
    required_namespaces = required_namespaces
  )
  class(wrapped_method) <- "dyneval::wrapped_method"
  wrapped_method
}

mt_objects <- list(
  create_method_type(
    name = "symmetric_distance_metric",
    input_types = list("x" = "matrix"),
    output_types = list("distance" = "symmetric_distance_matrix")
  ),
  create_method_type(
    name = "distance_metric",
    input_types = list("x" = "matrix", "y" = "matrix"),
    output_types = list("distance" = "distance_matrix")
  ),
  create_method_type(
    name = "symmetric_similarity_metric",
    input_types = list("x" = "matrix"),
    output_types = list("similarity" = "symmetric_similarity_matrix")
  ),
  create_method_type(
    name = "similarity_metric",
    input_types = list("x" = "matrix", "y" = "matrix"),
    output_types = list("similarity" = "similarity_matrix")
  ),
  create_method_type(
    name = "similarity_to_distance",
    input_types = list("similarity" = "similarity_matrix"),
    output_types = list("distance" = "distance_matrix")
  ),
  create_method_type(
    name = "distance_to_similarity",
    input_types = list("distance" = "distance_matrix"),
    output_types = list("similarity" = "similarity_matrix")
  ),
  create_method_type(
    name = "symmetric_similarity_to_distance",
    input_types = list("similarity" = "symmetric_similarity_matrix"),
    output_types = list("distance" = "symmetric_distance_matrix")
  ),
  create_method_type(
    name = "symmetric_distance_to_similarity",
    input_types = list("distance" = "symmetric_distance_matrix"),
    output_types = list("similarity" = "symmetric_similarity_matrix")
  ),
  create_method_type(
    name = "distance_dimensionality_reduction",
    input_types = list("distance" = "symmetric_distance_matrix"),
    output_types = list("space" = "space")
  ),
  create_method_type(
    name = "matrix_dimensionality_reduction",
    input_types = list("x" = "matrix"),
    output_types = list("space" = "space")
  ),
  create_method_type(
    name = "trajectory_inference",
    input_types = list("space" = "space"),
    output_types = list("pseudotime" = "pseudotime", "trajectory" = "trajectory")
  )
)
names(mt_objects) <- sapply(mt_objects, function(mto) mto$name)

is_method_type <- function(object) {
  object %in% names(mt_objects)
}

check_method_type <- function(object) {
  if (!is_method_type(object)) stop("method_type ", sQuote(object), " is unexpected")
}

#' Running a dyneval method
#'
#' @param wrapped_method The wrapped method
#' @param data The wrapped data to be passed
#' @param parameters Parameters of the method
#'
#' @export
run_method <- function(wrapped_method, data, parameters) {
  wrapped_method$method_function(data, parameters)
}
