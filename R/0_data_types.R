create_data_type <- function(name) {
  l <- list(name = name)
  class(l) <- "dyneval::data_type"
  l
}

inherits <- function(sub_data_type, ...) {
  super_data_types <- c(...)
  sub_data_type$super <- c(sub_data_type$super, super_data_types)
  sub_data_type
}

dt_objects <- list(
  create_data_type("matrix"),
  create_data_type("vector"),
  create_data_type("list"),
  create_data_type("symmetric"),
  create_data_type("expression") %>% inherits("matrix"),
  create_data_type("distance_matrix") %>% inherits("matrix"),
  create_data_type("similarity_matrix") %>% inherits("matrix"),
  create_data_type("symmetric_distance_matrix") %>% inherits("symmetric", "distance_matrix"),
  create_data_type("symmetric_similarity_matrix") %>% inherits("symmetric", "similarity_matrix"),
  create_data_type("reduced_space") %>% inherits(matrix)
)
names(dt_objects) <- sapply(dt_objects, function(dto) dto$name)

is_data_type <- function(object) {
  object %in% names(dt_objects)
}

check_data_type <- function(object) {
  if (!is_data_type(object)) stop("data_type ", sQuote(object), " is unexpected")
}

#' Wrap a data object
#'
#' @param data_types the data type object
#' @param data_object the data itself
#'
#' @export
wrap_data_object <- function(data_type, data_object) {
  check_data_type(data_type)
  l <- list(data_type = data_type, data_object = data_object)
  class(l) <- "dyneval::data_object"
  l
}

implements <- function(class, superclass) {
  check_data_type(class)
  check_data_type(superclass)
  if (class == superclass) {
    T
  } else {
    any(sapply(dt_objects[[class]]$super, implements, superclass))
  }
}

instanceof <- function(wrapped_data_object, superclass) {
  implements(wrapped_data_object$data_type, superclass)
}

#' Returns whether an object is a data_object
#'
#' @param object the object to test
#'
#' @export
is_wrapped_data_object <- function(object) {
  class(object) == "dyneval::data_object"
}

#' Unwrap a data_object
#'
#' @param wrapped_data_object a wrapped data_object
#'
#' @export
unwrap_data_object <- function(wrapped_data_object) {
  if (!is_wrapped_data_object(wrapped_data_object)) {
    stop("Object is not a wrapped data object")
  }
  wrapped_data_object$data_object
}
