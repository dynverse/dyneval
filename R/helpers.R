#' @export
to_tibble <- function(list_of_rows) {
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

  list_of_cols %>% as.tibble()
}

#' @export
load_datasets <- function() {
  datasets_info <- readRDS(paste0(.datasets_location, "/datasets.rds"))

  task_wrapped <- lapply(seq_len(nrow(datasets_info)), function(dataset_num) {
    dataset_id <- datasets_info$id[[dataset_num]]
    dataset <- dyngen::load_dataset(dataset_id, contents = dyngen::contents_dataset(experiment = F))

    with(dataset, dyneval::wrap_ti_task_data(
      ti_type = model$modulenetname,
      name = info$id,
      ids = rownames(counts),
      state_names = gs$milestone_names,
      state_net = gs$milestone_net,
      state_percentages = gs$milestone_percentages %>% slice(match(rownames(counts), id)) %>% gather(state, percentage, -id) %>% filter(percentage > 0), # temporary fix
      counts = counts
    ))
  })
  task_tib <- to_tibble(task_wrapped)
  task_tib %>% left_join(datasets_info, by = c("name"="id"))
}

# copied from assigninNamespace
my_assignin_namespace <- function (x, value, ns, pos = -1, envir = as.environment(pos)) {
  nf <- sys.nframe()
  if (missing(ns)) {
    nm <- attr(envir, "name", exact = TRUE)
    if (is.null(nm) || substr(nm, 1L, 8L) != "package:")
      stop("environment specified is not a package")
    ns <- asNamespace(substring(nm, 9L))
  }
  else ns <- asNamespace(ns)
  ns_name <- getNamespaceName(ns)
  if (bindingIsLocked(x, ns)) {
    in_load <- Sys.getenv("_R_NS_LOAD_")
    if (nzchar(in_load)) {
      if (in_load != ns_name) {
        msg <- gettextf("changing locked binding for %s in %s whilst loading %s",
                        sQuote(x), sQuote(ns_name), sQuote(in_load))
        if (!in_load %in% c("Matrix", "SparseM"))
          warning(msg, call. = FALSE, domain = NA, immediate. = TRUE)
      }
    }
    else if (nzchar(Sys.getenv("_R_WARN_ON_LOCKED_BINDINGS_"))) {
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
  }
  else {
    assign(x, value, envir = ns, inherits = FALSE)
  }
  if (!isBaseNamespace(ns)) {
    S3 <- .getNamespaceInfo(ns, "S3methods")
    if (!length(S3))
      return(invisible(NULL))
    S3names <- S3[, 3L]
    if (x %in% S3names) {
      i <- match(x, S3names)
      genfun <- get(S3[i, 1L], mode = "function", envir = parent.frame())
      if (.isMethodsDispatchOn() && methods::is(genfun,
                                                "genericFunction"))
        genfun <- methods::slot(genfun, "default")@methods$ANY
      defenv <- if (typeof(genfun) == "closure")
        environment(genfun)
      else .BaseNamespaceEnv
      S3Table <- get(".__S3MethodsTable__.", envir = defenv)
      remappedName <- paste(S3[i, 1L], S3[i, 2L], sep = ".")
      if (exists(remappedName, envir = S3Table, inherits = FALSE))
        assign(remappedName, value, S3Table)
    }
  }
  invisible(NULL)
}
