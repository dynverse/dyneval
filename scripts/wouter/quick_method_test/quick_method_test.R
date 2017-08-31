library(tidyverse)
library(dyneval)

tasks <- generate_toy_datasets()
tasks <- tasks[1, ]

counts <- tasks$counts[[1]]

# choose certain parameters for each method, at which we know this method will perform well for the toy dataset
method_descriptions <- list(
  waterfall=list(), # broken, cannot install RHmm on cluster due to lapack isue
  scorpius=list(),
  slingshot=list(),
  # # # slicer=list(max_same_milestone_distance=0.2, start_cell_id=progressions$percentage %>% which.min, min_branch_len=0.1, kmin=30, m=2), # broken it is inconceiveble really
  gpfates=list(nfates=1),
  stemid=list(clustnr=10, bootnr=10, pdishuf=10),
  tscan=list(),
  embeddr=list(nn_pct = 2),
  celltree_gibbs=list(sd_filter = 0),
  celltree_maptpx=list(sd_filter = 0),
  celltree_vem=list(sd_filter = 0),
  random_linear=list()
)

method_descriptions <- method_descriptions["gpfates"]

metric_names <- c("mean_R_nx", "auc_R_nx", "Q_local", "Q_global", "correlation", "ged", "isomorphic")

# test the methods and get the scores
results <- purrr::map(names(method_descriptions), function(method_name) {
# results <- parallel::mclapply(names(method_descriptions), function(method_name) {
# results <- PRISM::qsub_lapply(names(method_descriptions), function(method_name) {
  library(dyneval)

  cat("Processing ", method_name, "\n", sep="")
  method <- get(paste0("description_", method_name))()
  method_params <- method_descriptions[[method_name]]

  factory <- function(fun) {
    function(...) {
      warn <- err <- NULL
      res <- withCallingHandlers(
        tryCatch(fun(...), error=function(e) {
          err <<- conditionMessage(e)
          NULL
        }), warning=function(w) {
          warn <<- append(warn, conditionMessage(w))
          invokeRestart("muffleWarning")
        })
      list(res, warn=warn, err=err, method_name=method_name)
    }
  }

  factory(function() {
    dyneval:::execute_evaluation(tasks, method, method_params, metrics = metric_names, suppress_output = F)
  }
  )()
# }, mc.cores=8)
})
# }, qsub_environment = list2env(lst(method_descriptions, tasks, metric_names)), qsub_config = PRISM::override_qsub_config(memory = "10G"))

scores <- purrr::map(results, function(.) {
  bind_cols(
    tibble(error=list(.$err), warning=list(.$warning)),
    attr(.[[1]], "extras")$.summary,
    tibble(method_name=.$method_name)
  )
}) %>% bind_rows() %>% mutate(errored = !map_lgl(error, is.null))# %>% dplyr::select(-method_name1)

scores %>%
  dplyr::select(-starts_with("time")) %>%
  gather(score_id, score, -method_name, -method_short_name, -task_id, -error, -warning) %>%
  ggplot() + geom_bar(aes(method_name, score, fill=score_id), stat="identity") + facet_wrap(~score_id) + coord_flip()

