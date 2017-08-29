#' Attempts to convert a list of lists to a tibble
#'
#' @param list_of_rows The list to be converted to a tibble
#'
#' @return A tibble with the same number of rows as there were elements in \code{list_of_rows}
#' @export
#'
#' @examples
#' l <- list(list(a = 1, b = log10), list(a = 2, b = sqrt))
#' tib <- list_as_tibble(l)
#' tib
list_as_tibble <- function(list_of_rows) {
  list_names <- names(list_of_rows[[1]])

  list_of_cols <- lapply(seq_along(list_names), function(x) {
    colname <- list_names[[x]]
    list <- lapply(list_of_rows, function(z) z[[colname]])
    if (typeof(list[[1]]) != "list" && all(sapply(list, length) == 1)) {
      unlist(list, recursive = F)
    } else {
      list
    }
  }) %>% setNames(list_names)

  list_of_cols %>% as_tibble()
}

#' Extracts one row from a tibble and converts it to a list
#'
#' @param tib the tibble
#' @param row_id the index of the row to be selected
#'
#' @return the corresponding row from the tibble as a list
#' @export
#'
#' @examples
#' l <- list(list(a = 1, b = log10), list(a = 2, b = sqrt))
#' tib <- list_as_tibble(l)
#'
#' extract_row_to_list(tib, 2)
extract_row_to_list <- function(tib, row_id) {
  tib[row_id, ] %>% as.list %>% map(function(x) {
    if (is.null(x) | !is.list(x)) {
      x
    } else {
      x[[1]]
    }
  })
}


#' Assigning an object in a different namespace.
#'
#' It was copied from \code{\link[utils]{assignInNamespace}},
#' but a certain check was removed so that the set.seed function can be overwritten.
#'
#' This function is not to be used aside from overriding the set.seed function within the
#' TI evaluation.
my_assignin_namespace <- function (x, value, ns, pos = -1, envir = as.environment(pos)) {
  nf <- sys.nframe()
  if (missing(ns)) {
    nm <- attr(envir, "name", exact = TRUE)
    if (is.null(nm) || substr(nm, 1L, 8L) != "package:") {
      stop("environment specified is not a package")
    }
    ns <- asNamespace(substring(nm, 9L))
  } else {
    ns <- asNamespace(ns)
  }
  ns_name <- getNamespaceName(ns)
  # This check is removed!
  # if (nf > 1L) {
  #   if (ns_name %in% tools:::.get_standard_package_names()$base)
  #     stop("locked binding of ", sQuote(x), " cannot be changed",
  #          domain = NA)
  # }
  if (bindingIsLocked(x, ns)) {
    in_load <- Sys.getenv("_R_NS_LOAD_")
    if (nzchar(in_load)) {
      if (in_load != ns_name) {
        msg <- gettextf("changing locked binding for %s in %s whilst loading %s",
                        sQuote(x), sQuote(ns_name), sQuote(in_load))
        if (!in_load %in% c("Matrix", "SparseM"))
          warning(msg, call. = FALSE, domain = NA, immediate. = TRUE)
      }
    } else if (nzchar(Sys.getenv("_R_WARN_ON_LOCKED_BINDINGS_"))) {
      warning(gettextf("changing locked binding for %s in %s",
                       sQuote(x), sQuote(ns_name)), call. = FALSE,
              domain = NA, immediate. = TRUE)
    }
    unlockBinding(x, ns)
    assign(x, value, envir = ns, inherits = FALSE)
    w <- options("warn")
    on.exit(options(w))
    options(warn = -1)
    lockBinding(x, ns)
  } else {
    assign(x, value, envir = ns, inherits = FALSE)
  }
  if (!isBaseNamespace(ns)) {
    S3 <- .getNamespaceInfo(ns, "S3methods")
    if (!length(S3)) {
      return(invisible(NULL))
    }
    S3names <- S3[, 3L]
    if (x %in% S3names) {
      i <- match(x, S3names)
      genfun <- get(S3[i, 1L], mode = "function", envir = parent.frame())
      if (.isMethodsDispatchOn() && methods::is(genfun, "genericFunction")) {
        genfun <- methods::slot(genfun, "default")@methods$ANY
      }
      defenv <-
        if (typeof(genfun) == "closure") {
          environment(genfun)
        } else {
          .BaseNamespaceEnv
        }
      S3Table <- get(".__S3MethodsTable__.", envir = defenv)
      remappedName <- paste(S3[i, 1L], S3[i, 2L], sep = ".")
      if (exists(remappedName, envir = S3Table, inherits = FALSE))
        assign(remappedName, value, S3Table)
    }
  }
  invisible(NULL)
}
