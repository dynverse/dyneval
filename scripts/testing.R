library(tidyverse)
library(magrittr)
library(dyneval)

datasets <- readRDS("../dyngen/results/datasets.rds")
dataset <- datasets[[1]]
rm(datasets)

expr <- log2(dataset$counts+1)
wrapped_expr <- wrap_data_object(dyneval:::dt_matrix, expr)

# pick the first symmetric similarity metric and run it
wrapped_simmat <- run_method(dyneval:::imp_symmetric_similarity_metric[[1]], data = list(x = wrapped_expr), method = "pearson")

# pick the first sym2dist method and run it
wrapped_distmat <- run_method(dyneval:::imp_similarity_to_distance[[1]], data = wrapped_simmat, method = "corr")

# pick the first distance2space method and run it
wrapped_space <- run_method(dyneval:::imp_distance_to_space[[1]], data = wrapped_distmat, num_dimensions = 2)

# unwrap data and plot it
space <- unwrap_data_object(wrapped_space$space)
plot(space)



