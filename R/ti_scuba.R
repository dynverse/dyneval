#' @import ParamHelpers
#' @importFrom mclust mclust.options
#' @export
description_scuba <- function() {
  list(
    name = "SCUBA",
    short_name = "SCUBA",
    package_load = c("R.matlab", "readr"),
    package_installed = c(),
    par_set = makeParamSet(
      makeDiscreteParam(id = "cluster_mode", default = "pca2", values = c("original", "pca", "pca2"))
    ),
    properties = c(),
    run_fun = run_scuba,
    plot_fun = plot_scuba
  )
}

#' @export
run_scuba <- function(counts,
                      cluster_mode = "pca2") {
  # create new folder for data
  dataset_path <- tempfile()
  dir.create(dataset_path, recursive = T)

  # make a bunch of paths
  code_path <- paste0(path.package("dyneval"), "/extra_code/SCUBA")
  data_file <- paste0(dataset_path, "/Data.txt")
  projection_file <- paste0(dataset_path, "/intermediate_files/projection_all_data.txt")
  tree_file <- paste0(dataset_path, "/intermediate_files/final_tree.txt")
  bifurcation_file <- paste0(dataset_path, "/intermediate_files/bifurcation_direction.txt")
  treemat_file <- paste0(dataset_path, "/intermediate_files/final_tree.mat")

  # combine data
  combined_data <- data.frame(colnames(counts), t(counts), row.names = NULL)
  colnames(combined_data) <- c("Cell ID", rownames(counts))
  write.table(combined_data, data_file, sep = "\t", col.names = T, row.names = F)

  # run SCUBA
  command <- paste0(
    "cd '", code_path, "'; ",
    "module load matlab; ",
    "matlab -nodisplay -nodesktop -r \"", paste0(
      "addpath(genpath('", code_path, "/drtoolbox')); ",
      "RNAseq_preprocess('", dataset_path, "', 1, 1); ",
      "SCUBA('", dataset_path, "', '", cluster_mode, "'); ",
      "exit;"
    ), "\"")

  system(command)

  # read output
  projection_header <- readr::read_tsv(projection_file, n_max = 1, col_types = cols(Time = "c", .default = "d"))
  projection_data <- readr::read_tsv(projection_file, skip = 2, col_names = colnames(projection_header), col_types = cols(Time = "c", .default = "d"))
  tree_data <- readr::read_tsv(tree_file, col_types = cols(Time = "i", "Cluster ID" = "i", "Parent cluster" = "i", .default = "d"))
  bifurcation_header <- readr::read_tsv(bifurcation_file, n_max = 1, col_types = cols(Time = "c", .default = "d"))
  bifurcation_data <- readr::read_tsv(bifurcation_file, skip = 2, col_names = colnames(bifurcation_header), col_types = cols(Time = "c", .default = "d"))
  treemat_data <- R.matlab::readMat(treemat_file)$T[,,1]

  tree <- tree_data %>% select(cluster = `Cluster ID`, parent = `Parent cluster`, time_index = Time, cell_stage = `Cell Stage`) %>% mutate(cluster_name = paste0("Cluster", cluster))

  # create final output
  milestone_ids <- tree$cluster_name
  milestone_network <- tree %>%
    filter(parent != 0) %>%
    select(to = cluster_name, from_ix = parent) %>%
    left_join(tree %>% select(from = cluster_name, from_ix = cluster), by = "from_ix") %>%
    mutate(length = 1) %>%
    select(from, to, length)
  milestone_percentages <- data_frame(cell_id = projection_data$Time, milestone_id = paste0("Cluster", treemat_data$s[1,]), percentage = 1)

  # remove temporary output
  unlink(dataset_path, recursive = T)

  wrap_ti_prediction(
    ti_type = "linear",
    id = "SCUBA",
    cell_ids = rownames(counts),
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    milestone_percentages = milestone_percentages
    # ... add more output
  )
}

#' @export
plot_scuba <- function(ti_predictions) {
  stop("TODO")
  #qplot(percent_rank(ti_predictions$milestone_percentages[,1]), ti_predictions$milestone_percentages[,1], colour = data$sample_info$group.name)
}

