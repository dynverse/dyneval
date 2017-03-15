create_data_type <- function(data_types) {
  l <- list(data_types = data_types)
  class(l) <- "dyneval::data_type"
  l
}

inherits <- function(data_type1, data_type2) {
  combined_data_types <- unique(c(data_type1$data_types, data_type2$data_types))
  create_data_type(combined_data_types)
}

dt_matrix <- create_data_type("matrix")
dt_vector <- create_data_type("vector")
dt_list <- create_data_type("list")
dt_symmetric <- create_data_type("symmetric")

dt_distance_matrix <- create_data_type("distance_matrix") %>% inherits(dt_matrix)
dt_similarity_matrix <- create_data_type("similarity_matrix") %>% inherits(dt_matrix)
dt_symmetric_distance_matrix <- dt_symmetric %>% inherits(dt_distance_matrix)
dt_symmetric_similarity_matrix <- dt_symmetric %>% inherits(dt_similarity_matrix)

dt_reduced_space <- create_data_type("reduced_space") %>% inherits(dt_matrix)

wrap_data_object <- function(data_types, data_object) {
  l <- list(data_types = data_types$data_types, data_object = data_object)
  class(l) <- "dyneval::data_object"
  l
}

is_wrapped_data_object <- function(object) {
  class(object) == "dyneval::data_object"
}

unwrap_data_object <- function(wrapped_data_object) {
  if (!is_wrapped_data_object(wrapped_data_object)) {
    stop("Object is not a wrapped data object")
  }
  wrapped_data_object$data_object
}
