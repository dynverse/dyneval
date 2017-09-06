#' Description for SCUBA
#' @export
description_scuba <- function() create_description(
  name = "SCUBA",
  short_name = "SCUBA",
  package_loaded = c(),
  package_required = c("jsonlite", "readr"),
  par_set = makeParamSet(
    makeLogicalParam(id = "rigorous_gap_stats", default = T),
    makeIntegerParam(id = "N_dim", lower=2, upper=20, default=2),
    makeNumericParam(id = "low_gene_threshold", lower = 0, upper = 5, default = 1),
    makeNumericParam(id = "low_gene_fraction_max", lower = 0, upper = 1, default = 0.7),
    makeIntegerParam(id = "min_split", lower=1, upper=100, default=15),
    makeNumericParam(id = "min_percentage_split", lower=0, upper=1, default=0.25)
  ),
  properties = c(),
  run_fun = run_scuba,
  plot_fun = plot_scuba,
  make_command = paste0("extra_code/PySCUBA/make ", get_dyneval_install_path(), "/scuba")
)

#' @importFrom utils write.table
run_scuba <- function(counts,
                      rigorous_gap_stats=TRUE,
                      N_dim=2,
                      low_gene_threshold=1,
                      low_gene_fraction_max=0.7,
                      min_split=15,
                      min_percentage_split=0.25) {

  requireNamespace("jsonlite")
  temp_folder <- tempfile()
  dir.create(temp_folder, recursive = T)

  counts %>%
    #{log2(. + 1)} %>%
    t %>%
    as.data.frame() %>%
    write.table(paste0(temp_folder, "/counts.tsv"), sep="\t", col.names=NA)

  system2(
    "/bin/bash",
    args = c(
      "-c",
      shQuote(glue::glue(
        "cd {get_scuba_path()}/",
        "source bin/activate",
        "python3 {find.package('dyneval')}/extra_code/PySCUBA/wrapper.py {temp_folder} 0 {c(0, 1)[as.numeric(rigorous_gap_stats)+1]} {N_dim} {low_gene_threshold} {low_gene_fraction_max} {min_split} {min_percentage_split}",
        .sep = ";"))
    )
  )

  output <- jsonlite::read_json(paste0(temp_folder, "/output.json"))
  labels <- as.character(unlist(output$labels))

  milestone_network <- output$new_tree %>%
    map(unlist) %>%
    .[-1] %>%
    map(as.numeric) %>%
    do.call(rbind, .) %>%
    .[, c(1, 3)] %>%
    as.data.frame() %>%
    magrittr::set_colnames(c("from", "to")) %>%
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
