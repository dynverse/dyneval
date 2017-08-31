library(cowplot)
library(tidyverse)
library(dyneval)
library(mlrMBO)
library(parallelMap)

tasks <- generate_toy_datasets()
counts <- tasks$counts[[1]]

## choose a method
# method <- description_scuba()
# method <- description_waterfall()
# method <- description_monocle_ddrtree()
# method <- description_scorpius()
method <- description_slingshot()
# method <- description_scuba()
method_fun <- make_obj_fun(method)

## apply default params
default_params <- generateDesignOfDefaults(method$par_set)
params <- generateDesign(n = 20, par.set = method$par_set)

#eval_out <- method_fun(default_params[1,] %>% as.list, tasks = tasks)
eval_out <- method_fun(list(), tasks = tasks[1:min(nrow(tasks), 10),])
eval_extras <- attr(eval_out,"extras")
attr(eval_out,"extras") <- NULL
eval_out
eval_extras$.summary

## apply it manually
method_out <- method$run_fun(counts)

corank_out <- compute_coranking(tasks$geodesic_dist[[1]], method_out$geodesic_dist)





model <- eval_extras$.models[[1]]
task <- tasks %>% extract_row_to_list(1)

milpct_1 <- model$milestone_percentages %>% reshape2::acast(cell_id ~ milestone_id, value.var = "percentage", fill = 0, drop = F)
milpct_2 <- task$milestone_percentages %>% reshape2::acast(cell_id ~ milestone_id, value.var = "percentage", fill = 0, drop = F)
testthat::expect_equal(rownames(milpct_1), rownames(milpct_2))

node_cor <- (cor(milpct_1, milpct_2) + 1) / 2

# try GraphAlignment
library(GraphAlignment)
pinitial <- InitialAlignment(psize = 3, r = node_cor, mode="reciprocal")

lookupLink <- seq(0, 1, .5)

adj_1 <- model$milestone_network %>%
  mutate(from = factor(from, levels = model$milestone_ids), to = factor(to, levels = model$milestone_ids)) %>%
  reshape2::acast(from ~ to, value.var = "length", fill = Inf, drop = F)
adj_1 <- pmin(adj_1, t(adj_1))
adj_2 <- task$milestone_network %>%
  mutate(from = factor(from, levels = task$milestone_ids), to = factor(to, levels = task$milestone_ids)) %>%
  reshape2::acast(from ~ to, value.var = "length", fill = Inf, drop = F)
adj_2 <- pmin(adj_2, t(adj_2))

linkParams <- ComputeLinkParameters(1 / adj_1, 1 / adj_2, pinitial, lookupLink)

lookupNode<-c(-.5,.5,1.5)
nodeParams <- ComputeNodeParameters(
  dimA = length(model$milestone_ids),
  dimB = length(task$milestone_ids),
  node_cor,
  pinitial,
  lookupNode)

al<-AlignNetworks(A=adj_1, B=adj_2, R=node_cor, P=pinitial,
                  linkScore=linkParams$ls,
                  selfLinkScore=linkParams$ls,
                  nodeScore1=nodeParams$s1, nodeScore0=nodeParams$s0,
                  lookupLink=lookupLink, lookupNode=lookupNode,
                  bStart=.1, bEnd=30,
                  maxNumSteps=50)

al
#
# # sudo apt-get install libglpk-dev
# library(Corbi)
#
#
# # library(GraphAlignment)
# # pinitial <- InitialAlignment(psize=34, r=ex$r, mode="reciprocal")


library(netdist)
?netdist::net_emd
