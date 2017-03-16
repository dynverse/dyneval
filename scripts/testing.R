library(tidyverse)
library(magrittr)
library(dyneval)

datasets <- readRDS("../dyngen/results/datasets.rds")
dataset <- datasets[[1]]
rm(datasets)

expr <- log2(dataset$counts+1)
wrapped_expr <- wrap_data_object("expression", expr)
named_wrapped_expr <- list(x = wrapped_expr)

# pick the first symmetric similarity metric and run it
meth1 <- dyneval:::imp_symmetric_similarity_metric[[1]]
wrapped_simmat <- run_method(meth1, named_wrapped_expr, get_parameter_row(meth1, 1))

# pick the first sym2dist method and run it
meth2 <- dyneval:::imp_symmetric_similarity_to_distance[[1]]
wrapped_distmat <- run_method(meth2, wrapped_simmat, get_parameter_row(meth2, 1))

# pick the first distance2space method and run it
meth3 <- dyneval:::imp_distance_to_space[[1]]
wrapped_space <- run_method(meth3, wrapped_distmat, get_parameter_row(meth3, 1))

meth4 <- dyneval:::imp_trajectory_inference[[1]]
wrapped_traj <- run_method(meth4, wrapped_space, get_parameter_row(meth4, 1))

# unwrap data and plot it
space <- unwrap_data_object(wrapped_space$space)
plot(space)

## making a method graph ------------------------------------------------------
mt_object_names <- names(dyneval:::mt_objects)
dt_object_names <- names(dyneval:::dt_objects)
mt_objects <- dyneval:::mt_objects
dt_objects <- dyneval:::dt_objects
mt_as_df <- bind_rows(lapply(mt_objects, function(mt_obj) {
  bind_rows(
    data_frame(name = mt_obj$name, param_name = names(mt_obj$input_types), param_type = unlist(mt_obj$input_types), type = "input"),
    data_frame(name = mt_obj$name, param_name = names(mt_obj$output_types), param_type = unlist(mt_obj$output_types), type = "output")
  )
}))
inherits_as_df <- bind_rows(lapply(dt_objects, function(dt_obj) {
  if (!is.null(dt_obj$super)) {
    data_frame(name = dt_obj$name, super = dt_obj$super, type = "inherits")
  } else {
    NULL
  }
}))
mt_as_df
object_types <- bind_rows(
  data_frame(id = mt_object_names, node_type = "method"),
  data_frame(id = dt_object_names, node_type = "data")
)

edge_df <- bind_rows(
  mt_as_df %>%
    mutate(from = ifelse(type == "input", param_type, name), to = ifelse(type == "input", name, param_type)) %>%
    select(from, to, label = param_name, name, param_name, param_type, type),
  inherits_as_df %>% mutate(from = name, to = super, label = type) %>% select(from, to, label, name, super, type)
)
type_colours <- c(input = "blue", output = "red", inherits = "lightgray")
node_df <- object_types %>% mutate(vertex.shape = ifelse(node_type == "method", "rectangle", "circle")) %>% select(id, vertex.shape)
library(igraph)
igr <- igraph::graph_from_data_frame(edge_df %>% as.data.frame(), directed = T, vertices = node_df)

pdf("scratch/methods.pdf", 10, 10)
plot(igr, vertex.shape = node_df$vertex.shape, edge.color = type_colours[edge_df$type])
dev.off()
write_tsv(edge_df, "scratch/methods_edges.tsv")
write_tsv(node_df, "scratch/methods_nodes.tsv")

## toying with params ---------------------------------------------------------
wrapped_method <- dyneval:::imp_similarity_to_distance[[1]]
# wrapped_method$parameter_sets <- list(
#   list(method = "alpha_method", alpha = c(.1, .2, .5), k = 2:4, fun = function(a) a + 1),
#   list(method = "beta_method", beta = c("one", "two"), k = 3:5, fun = function(b) b + 2)
# )
wrapped_method$parameter_sets <- list(
  list(method = "alpha_method", alpha = c(.1, .2, .5), k = 2:4),
  list(method = "beta_method", beta = c("one", "two"), k = 3:5)
)

## toying with mlr paramsets --------------------------------------------------
library(mlr)
params <- makeParamSet(
  makeIntegerLearnerParam(id = "ntree", default = 500L, lower = 1L),
  makeIntegerLearnerParam(id = "se.ntree", default = 100L, lower = 1L, when = "both", requires = quote(se.method == "bootstrap")),
  makeDiscreteLearnerParam(id = "se.method", default = "jackknife",
                           values = c("bootstrap", "jackknife",  "sd"),
                           requires = quote(se.method %in% c("jackknife") && keep.inbag == TRUE),
                           when = "both"),
  makeIntegerLearnerParam(id = "se.boot", default = 50L, lower = 1L, when = "both"),
  makeIntegerLearnerParam(id = "mtry", lower = 1L),
  makeLogicalLearnerParam(id = "replace", default = TRUE),
  makeUntypedLearnerParam(id = "strata", tunable = FALSE),
  makeIntegerVectorLearnerParam(id = "sampsize", lower = 1L),
  makeIntegerLearnerParam(id = "nodesize", default = 5L, lower = 1L),
  makeIntegerLearnerParam(id = "maxnodes", lower = 1L),
  makeLogicalLearnerParam(id = "importance", default = FALSE),
  makeLogicalLearnerParam(id = "localImp", default = FALSE),
  makeIntegerLearnerParam(id = "nPerm", default = 1L),
  makeLogicalLearnerParam(id = "proximity", default = FALSE, tunable = FALSE),
  makeLogicalLearnerParam(id = "oob.prox", requires = quote(proximity == TRUE), tunable = FALSE),
  makeLogicalLearnerParam(id = "do.trace", default = FALSE, tunable = FALSE),
  makeLogicalLearnerParam(id = "keep.forest", default = TRUE, tunable = FALSE),
  makeLogicalLearnerParam(id = "keep.inbag", default = FALSE, tunable = FALSE)
)

# work with tuneParms

# saving benchmarking results in data base? http://stackoverflow.com/questions/3285307/save-r-plot-to-web-server
