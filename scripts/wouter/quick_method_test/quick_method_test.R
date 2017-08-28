library(tidyverse)
library(dyneval)

source("scripts/wouter/toy/generation.R")
source("scripts/wouter/toy/perturbation.R")

# generate a simple linear toy dataset
task <- generate_linear()

timepoints <- task$progressions$percentage
cell_ids <- task$progressions$cell_id

ngenes <- 200
coefs <- runif(ngenes, -1, 1)
intercepts <- runif(ngenes, -2, 2)

expression <- t(coefs %*% t(timepoints)) %>% apply(1, function(x) x + intercepts) %>% t
expression <- expression + rnorm(length(expression), mean = 0, sd=0.1)
dimnames(expression) <- list(cell_ids, seq_len(ncol(expression)))
expression[order(timepoints),] %>% t %>% pheatmap::pheatmap(cluster_cols=F, scale="row")
plot(timepoints, expression[, 1])

counts <- round((2^expression) * 10)
task$counts <- counts
tasks <- dyneval:::list_as_tibble(list(task))

# choose certain parameters for each method, at which we know this method will perform well for the toy dataset
method_descriptions <- list(
  #random_linear=list(),
  # waterfall=list(), # broken
  scorpius=list(),
  slingshot=list(),
  #  slicer=list(max_same_milestone_distance=0.2, start_cell_id=progressions$percentage %>% which.min, min_branch_len=0.1, kmin=30, m=2), # broken it is inconceiveble really
  # gpfates=list(nfates=1),
  stemid=list(),
  tscan=list(),
  embeddr=list(),
  celltree_gibbs=list(sd_filter = 0),
  celltree_maptpx=list(sd_filter = 0),
  celltree_vem=list(sd_filter = 0)
)

# test the methods and get the scores
metric_names <- c("mean_R_nx", "auc_R_nx", "Q_local", "Q_global", "correlation", "ged", "isomorphic")

scores <- purrr::map(names(method_descriptions), function(method_name) {
  method <- get(paste0("description_", method_name))()
  method_params <- method_descriptions[[method_name]]

  walk(method$package_load, ~require(., character.only=TRUE))

  method_out <- dyneval:::execute_evaluation(tasks, method, method_params, metrics = metric_names)
  out_extras <- attr(method_out, "extras")

  out_extras$.summary %>% mutate(method_name = method_name)
}) %>% bind_rows()

scores$correlation

scores %>%
  dplyr::select(-starts_with("time")) %>%
  gather(score_id, score, -method_name) %>%
  ggplot() + geom_bar(aes(method_name, score), stat="identity") + facet_wrap(~score_id)
