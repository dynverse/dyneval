#' @export
install_gpfates <- function() {
  system(glue::glue(
    "bash ",
    path.package("dyneval"), "/extra_code/GPfates/make ",path.package("dyneval"), "/extra_code/GPfates/"
  ))
}

#' @import ParamHelpers
#' @export
description_gpfates <- function() {
  if(!dir.exists(glue::glue("{path.package('dyneval')}/extra_code/GPfates/gpfates"))) {
    warning("gpfates not installed, installing now")
    install_gpfates()
  }

  list(
    name = "GPfates",
    short_name = "GPfates",
    package_installed = c(),
    par_set = makeParamSet(
      makeNumericParam(id = "log_expression_cutoff", lower = 0.5, upper = 5, default = 2),
      makeNumericParam(id = "min_cells_expression_cutoff", lower = 0, upper = 20, default = 2),
      makeIntegerParam(id="nfates", lower=1, upper=20, default=1),
      makeIntegerParam(id="ndims", lower=1, upper=5, default=2)
    ),
    properties = c(),
    run_fun = run_gpfates,
    plot_fun = plot_gpfates
  )
}

# library(tidyverse)
# library(dyneval)
#
# .version <- 4
# .datasets_location <- paste0("/home/wouters/thesis/projects/dyngen/results/", .version, "/")
# dataset <- load_datasets(ndatasets = 1) %>% dyneval:::extract_row_to_list(1)
#
# counts <- dataset$counts


## TODO: give simulationtime as prior
run_gpfates <- function(
  counts,
  log_expression_cutoff=2,
  min_cells_expression_cutoff=2,
  nfates=2,
  ndims=2
) {

  temp_folder <- tempdir()

  counts %>%
    t %>%
    write.table(paste0(temp_folder, "expression.csv"), sep="\t")

  counts %>%
    {data.frame(cell_id=rownames(.), row.names=rownames(.))} %>%
    write.table(paste0(temp_folder, "cellinfo.csv"), sep="\t")

  system(glue::glue(
    "cd {path.package('dyneval')}/extra_code/GPfates/gpfates",
    "source bin/activate",
    "python3 ..//gpfates_wrapper.py {temp_folder} {log_expression_cutoff} {min_cells_expression_cutoff} {nfates} {ndims}"
,
  .sep = ";"))

  pseudotime <- read_csv(glue::glue("{temp_folder}pseudotimes.csv"), col_names = c("cell_id", "time"))

  phi <- read_csv(glue::glue("{temp_folder}phi.csv"), col_names = c("cell_id", glue::glue("M{seq_len(nfates)}")), skip = 1)

  dr <- read_csv(glue::glue("{temp_folder}dr.csv"), col_names = c("cell_id", glue::glue("Comp{seq_len(5)}")), skip = 1)

  pseudotime <- pseudotime %>% mutate(time=(time-min(time))/(max(time) - min(time)))

  # first get percentages of all final milestones, by getting their phi values, and multiplying by the pseudotime
  milestone_percentages <- phi %>%
    gather(milestone_id, percentage, -cell_id) %>%
    left_join(pseudotime, by="cell_id") %>%
    mutate(percentage=percentage*time)

  # now add the starting milestone, just 1-pseudotime
  milestone_percentages <- milestone_percentages %>% bind_rows(pseudotime %>% mutate(milestone_id="M0", percentage= 1-time)) %>% select(-time)

  # now create the other objects
  milestone_network <- tibble(
    from="M0",
    to=glue::glue("M{seq_len(nfates)}"),
    length=1
  )
  milestone_ids <- unique(c(milestone_network$from, milestone_network$to))

  wrap_ti_prediction(
    ti_type = "GPfates",
    id = "GPfates",
    cell_ids = rownames(counts),
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    milestone_percentages = milestone_percentages,
    dimred_samples = dr,
    pseudotime = pseudotime
  )
}

## TODO extract OMGP
#' @export
plot_gpfates <- function(ti_predictions) {
  sample_df <- data.frame(
    ti_predictions$dimred_samples
  ) %>% left_join(pseudotime, by="cell_id")
  ggplot() +
    geom_point(aes(Comp1, Comp2, colour = time), sample_df) +
    coord_equal() +
    viridis::scale_color_viridis()
}
