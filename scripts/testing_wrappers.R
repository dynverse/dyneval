library(dyneval)
library(dynutils)
library(dplyr)
library(tidyr)
library(purrr)
library(magrittr)
library(tibble)
library(ggplot2)

dataset <- dynutils::extract_row_to_list(dyntoy::toy_tasks, 5)
counts <- dataset$counts
# start_cell_id <- dataset$special_cells$start_cell_id
# cell_grouping <- dataset$cell_grouping

# dpt
start_cell_id = NULL
sigma = "local"
distance = "euclidean"
n_eigs = 20
density_norm = TRUE
n_local_lower = 5
n_local_upper = 7
w_width = .1

# embeddr
kernel = "nn"
metric = "correlation"
nn_pct = 1
eps = 1
t = 1
symmetrize = "mean"
measure_type = "unorm"
p = 2
thresh = .001
maxit = 10
stretch = 2
smoother = "smooth.spline"

# gpfates
nfates = 1
ndims = 2
log_expression_cutoff = 2
min_cells_expression_cutoff = 2
num_cores = 1
verbose = FALSE

# mfa
b = 2
iter=2000
thin=1
zero_inflation=FALSE
pc_initialise=1
prop_collapse=0
scale_input=TRUE

# monocle
reduction_method <- "DDRTree"
max_components <- 2
norm_method <- "vstExprs"
auto_param_selection <- TRUE
num_paths <- NULL

# monocle 1
reduction_method <- "ICA"
max_components <- 2
norm_method <- "vstExprs"
auto_param_selection <- TRUE
num_paths <- 1

# mpath
distMethod = "euclidean"
method = "kmeans"
numcluster = 11
diversity_cut = .6
size_cut = .05

# ouija
iter = 20
response_type = "switch"
inference_type = "hmc"
normalise_expression = TRUE

# phenopath
thin = 40
z_init = 1
model_mu = FALSE
scale_y = TRUE


