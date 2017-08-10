#' Converts a list of lists to a tibblee.
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

  list_of_cols %>% as.tibble()
}

#' @export
load_datasets <- function(mc_cores = 1) {
  datasets_info <- readRDS(paste0(.datasets_location, "/datasets.rds"))

  task_wrapped <- parallel::mclapply(seq_len(nrow(datasets_info)), mc.cores = mc_cores, function(dataset_num) {
    dataset_id <- datasets_info$id[[dataset_num]]
    dataset <- dyngen::load_dataset(dataset_id)

    list2env(dataset, environment())

    dyneval::wrap_ti_task_data(
      ti_type = model$modulenetname,
      id = dataset_id,
      cell_ids = rownames(counts),
      milestone_ids = gs$milestone_percentages$milestone %>% unique %>% as.character, # change this to milestone_id in the future?
      milestone_network = gs$milestone_network %>%
        mutate(
          from = as.character(from),
          to = as.character(to),
          length = ifelse(!is.na(length), length, 1)
        ),
      milestone_percentages = cellinfo %>%
        left_join(gs$milestone_percentages, by = c("step_id"="cell_id")) %>%
        filter(percentage > 0) %>%
        rename(milestone_id = milestone) %>% # change this to milestone_id
        mutate(milestone_id = as.character(milestone_id)) %>%
        select(cell_id, milestone_id, percentage),
      counts = counts,
      sample_info = cellinfo %>%
        select(cell_id, step, simulation_time = simulationtime),
      platform_id = platform$platform_id,
      experiment_type = gsub(".*_[0-9]+_([^_]*)_.*", "\\1", dataset_id) # temporary fix
    )
  })
  task_wrapped %>%
    list_as_tibble %>%
    left_join(datasets_info, by = c("id" = "id")) %>%
    mutate(task_ix = seq_len(n())) %>%
    group_by(ti_type) %>%
    mutate(subtask_ix = seq_len(n())) %>%
    ungroup
}

#' This function assigns an object in a different namespace.
#'
#' It was copied from \code{\link[utils]{assignInNamespace}},
#' but without an extra check so that I can override the set.seed function.
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
