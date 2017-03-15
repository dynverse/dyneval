create_method_type <- function(name, input_types, output_types) {
  l <- list(
    name = name,
    input_types = input_types,
    output_types = output_types
  )
  class(l) <- "dyneval::method_type"
  l
}

wrap_method <- function(method_name, method_type, method_function, parameter_sets, required_namespaces) {
  input_types <- method_type$input_types
  output_types <- method_type$output_types
  form_args <- formalArgs(method_function)
  if (!all(names(input_types) == form_args[seq_along(input_types)])) {
    stop("Not all input_types were found")
  }
  wrapped_function <- function(...) {
    execution_arguments <- list(...)
    for (namespace in required_namespaces) {
      requireNamespace(namespace)
    }
    for (i in seq_along(input_types)) {
      if (!all(input_types[[i]]$data_types %in% execution_arguments[[i]]$data_types)) {
        message <- paste0(
          "Mismatch in input data types.\n",
          "Parameter: ", names(execution_arguments)[[i]], "\n",
          "expected DT: ", paste(input_types[[i]]$data_types, collapse=";"), "\n",
          "actual DT: ", paste(execution_arguments[[i]]$data_types, collapse=";"))
        stop(message)
      }
      execution_arguments[[i]] <- unwrap_data_object(execution_arguments[[i]])
    }
    method_output <- do.call(method_function, execution_arguments)
    if (!all(names(output_types) == names(method_output))) {
      stop("Not all output_types were found")
    }
    for (i in seq_along(output_types)) {
      if (!all(output_types[[i]]$data_types %in% method_output[[i]]$data_types)) {
        message <- paste0(
          "Mismatch in output data types.\n",
          "Parameter: ", names(method_output)[[i]], "\n",
          "expected DT: ", paste(output_types[[i]]$data_types, collapse=";"), "\n",
          "actual DT: ", paste(method_output[[i]]$data_types, collapse=";"))
        stop(message)
      }
    }
    method_output
  }
  wrapped_method <- list(
    method_type_name = method_type$name,
    method_name = method_name,
    method_function = wrapped_function,
    parameter_sets = parameter_sets,
    required_namespaces = required_namespaces
  )
  class(wrapped_method) <- "dyneval::wrapped_method"
  wrapped_method
}

mt_symmetric_distance_metric <- create_method_type(
  name = "symmetric_distance_metric",
  input_types = list("x" = dt_matrix),
  output_types = list("distance" = dt_symmetric_distance_matrix)
)
mt_distance_metric <- create_method_type(
  name = "nonsymmetric_distance_metric",
  input_types = list("x" = dt_matrix, "y" = dt_matrix),
  output_types = list("distance" = dt_distance_matrix)
)
mt_symmetric_similarity_metric <- create_method_type(
  name = "symmetric_similarity_metric",
  input_types = list("x" = dt_matrix),
  output_types = list("similarity" = dt_symmetric_similarity_matrix)
)
mt_similarity_metric <- create_method_type(
  name = "nonsymmetric_similarity_metric",
  input_types = list("x" = dt_matrix, "y" = dt_matrix),
  output_types = list("similarity" = dt_similarity_matrix)
)
mt_similarity_to_distance <- create_method_type(
  name = "similarity_to_distance",
  input_types = list("similarity" = dt_similarity_matrix),
  output_types = list("distance" = dt_distance_matrix)
)
mt_distance_to_similarity <- create_method_type(
  name = "distance_to_similarity",
  input_types = list("distance" = dt_distance_matrix),
  output_types = list("similarity" = dt_similarity_matrix)
)
mt_distance_dimensionality_reduction <- create_method_type(
  name = "distance_dimensionality_reduction",
  input_types = list("distance" = dt_distance_matrix),
  output_types = list("space" = dt_reduced_space)
)
mt_matrix_dimensionality_reduction <- create_method_type(
  name = "matrix_dimensionality_reduction",
  input_types = list("x" = dt_matrix),
  output_types = list("space" = dt_reduced_space)
)

#' Running a dyneval method
#'
#' @param wrapped_method The wrapped method
#' @param data The wrapped data to be passed
#' @param ... Parameters of the method
#'
#' @export
run_method <- function(wrapped_method, data, ...) {
  params <- list(...)
  combined_data <- c(
    data,
    params
  )
  do.call(wrapped_method$method_function, combined_data)
}
