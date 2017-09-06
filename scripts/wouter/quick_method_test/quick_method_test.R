library(tidyverse)
library(dyneval)

set.seed(1)
tasks <- generate_toy_datasets()
tasks <- tasks[1,]

counts <- tasks$counts[[1]]
special_cells <- tasks$special_cells[[1]]
cell_grouping <- tasks$cell_grouping[[1]]
progressions <- tasks$progressions[[1]]


# choose certain parameters for each method, at which we know this method will perform well for the toy dataset
method_descriptions <- list(
  waterfall=list(),
  scorpius=list(),
  slingshot=list(),
  gpfates=list(nfates=1),
  stemid=list(clustnr=10, bootnr=10, pdishuf=10),
  tscan=list(),
  embeddr=list(nn_pct = 2),
  celltree_gibbs=list(sd_filter = 0),
  celltree_maptpx=list(sd_filter = 0),
  celltree_vem=list(sd_filter = 0),
  scuba=list(),
  slicer = list(min_branch_len=50),
  monocle_ddrtree=list(),
  wishbone = list(branch=F),
  pseudogp = list(iter=100, initialise_from="principal_curve"),
  mpath = list(numcluster=6),
  random_linear=list()
)
#
#
# prediction <- dyneval:::execute_method(tasks, description_mpath(), method_descriptions$mpath, suppress_output = F)[[1]]$model

# method_descriptions <- method_descriptions["mpath"]

metric_names <- c("mean_R_nx", "auc_R_nx", "Q_local", "Q_global", "correlation", "isomorphic", "robbie_network_score")

# test the methods and get the scores
results <- purrr::map(names(method_descriptions), function(method_name) {
# results <- pbapply::pblapply(names(method_descriptions), cl = 7, function(method_name) {
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
# }, qsub_environment = list2env(lst(method_descriptions, tasks, metric_names)), qsub_config = PRISM::override_qsub_config(memory = "10G", execute_before=c("module load python/x86_64/3.5.1")))
#
# PRISM::qsub_run(
#   function(i) {dyneval:::run_slicer(counts)},
#   qsub_environment=list2env(lst(counts)),
#   qsub_config = PRISM::override_qsub_config(
#     memory = "10G",
#     execute_before=c("module load python/x86_64/3.5.1"),
#     r_module = "R",
#     remove_tmp_folder=F
#   )
# )

scores <- purrr::map(results, function(.) {
  bind_cols(
    tibble(error=list(.$err), warning=list(.$warning)),
    attr(.[[1]], "extras")$.summary,
    tibble(method_name=.$method_name)
  )
}) %>% bind_rows() %>% mutate(errored = !map_lgl(error, is.null)) %>% dplyr::select(-method_name1)

scores %>%
  dplyr::select(-starts_with("time")) %>%
  gather(score_id, score, -method_name, -method_short_name, -task_id, -error, -warning) %>%
  ggplot() + geom_bar(aes(method_name, score, fill=score_id), stat="identity") + facet_wrap(~score_id) + coord_flip()

