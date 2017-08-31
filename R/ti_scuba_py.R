get_scuba_path <- function() {
  "~/.dyneval/scuba/"
}

#' Installing SCUBA
#' @export
install_scuba <- function() {
  dir.create(get_scuba_path(), showWarnings = FALSE, recursive = TRUE)

  system(glue::glue(
    "bash ",
    path.package("dyneval"), "/extra_code/PySCUBA/make ", get_scuba_path()
  ))
}

#' Description for SCUBA
#' @export
description_scuba <- function() {
  list(
    name = "SCUBA",
    short_name = "SCUBA",
    package_load = c(),
    package_installed = c("jsonlite", "readr"),
    par_set = makeParamSet(
      makeLogicalParam(id = "log_mode"),
      makeLogicalParam(id = "rigorous_gap_stats"),
      makeIntegerParam(id = "N_dim", lower=2, upper=20, default=2),
      makeNumericParam(id = "low_gene_threshold", lower = 0, upper = 5, default = 1),
      makeNumericParam(id = "low_gene_fraction_max", lower = 0, upper = 1, default = 0.7),
      makeIntegerParam(id = "min_split", lower=1, upper=100, default=15),
      makeNumericParam(id = "min_percentage_split", lower=0, upper=1, default=0.25),
    ),
    properties = c(),
    run_fun = run_scuba,
    plot_fun = plot_scuba
  )
}

#' @export
#' @importFrom utils write.table
run_scuba <- function(counts,
                      log_mode=TRUE,
                      rigorous_gap_stats=TRUE,
                      N_dim=2,
                      low_gene_threshold=1,
                      low_gene_fraction_max=0.7,
                      min_split=15,
                      min_percentage_split=0.25) {

  temp_folder <- tempfile()
  dir.create(temp_folder, recursive = T)

  counts %>%
    t %>%
    as.data.frame() %>%
    write.table(paste0(temp_folder, "/counts.tsv"), sep="\t", col.names=NA)

  system2(
    "/bin/bash",
    args = c(
      "-c",
      shQuote(glue::glue(
        "cd {get_scuba_path()}/scuba/",
        "source bin/activate",
        "python {path.package('dyneval')}/extra_code/PySCUBA/scuba_wrapper.py {temp_folder} {c(0, 1)[as.numeric(log_mode)]} {c(0, 1)[as.numeric(rigorous_gap_stats)]} {N_dim} {low_gene_threshold} {low_gene_fraction_max} {min_split} {min_percentage_split}",
        .sep = ";"))
    )
  )

  output <- jsonlite::read_json(paste0(temp_folder, "/output.json"))
  labels <- as.character(unlist(output$labels))
  milestone_network <- readr::read_tsv(paste0(temp_folder, "/final_tree.tsv")) %>%  # tree was outputted seperatly by a python function, this is not returned by the function but can only be extracted by saving it
    rename(to=`Cluster ID`, stage=Stage, from=`Parent cluster`) %>%
    select(from, to) %>%
    mutate(
      length=1,
      from=as.character(from),
      to=as.character(to)
    )

  milestone_ids <- unique(c(milestone_network$from, milestone_network$to))

  milestone_percentages <- data_frame(
    cell_id = rownames(counts),
    milestone_id = labels,
    percentage = 1
  )

  # remove temporary output
  unlink(temp_folder, recursive = T)

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

plot_scuba <- function(ti_predictions) {
  stop("TODO")
  #qplot(percent_rank(ti_predictions$milestone_percentages[,1]), ti_predictions$milestone_percentages[,1], colour = data$sample_info$group.name)
}
