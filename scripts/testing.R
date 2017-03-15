library(tidyverse)
library(magrittr)
library(dyneval)

datasets <- readRDS("../dyngen/results/datasets.rds")
dataset <- datasets[[1]]
rm(datasets)

expr <- log2(dataset$counts+1)
wrapped_expr <- dyneval:::wrap_data_object(dyneval:::dt_matrix, expr)

wrapped_simmat <- dyneval:::imp_symmetric_similarity_metric[[1]]$method_function(x = wrapped_expr, method = "pearson")
wrapped_distmat <- dyneval:::imp_similarity_to_distance[[1]]$method_function(similarity = wrapped_simmat$similarity, method = "corr")
wrapped_space <- dyneval:::imp_distance_to_space[[1]]$method_function(distance = wrapped_distmat$distance, num_dimensions = 2)
space <- dyneval:::unwrap_data_object(wrapped_space$space)

plot(space)
