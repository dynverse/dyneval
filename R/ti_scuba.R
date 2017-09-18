#' Description for SCUBA
#' @export
description_scuba <- function() create_description(
  name = "SCUBA",
  short_name = "SCUBA",
  package_loaded = c(),
  package_required = c("jsonlite", "readr", "SCUBA"),
  par_set = makeParamSet(
    makeLogicalParam(id = "rigorous_gap_stats", default = TRUE),
    makeIntegerParam(id = "N_dim", lower=2L, upper=3L, default=2L),
    makeNumericParam(id = "low_gene_threshold", lower = 0, upper = 5, default = 1),
    makeNumericParam(id = "low_gene_fraction_max", lower = 0, upper = 1, default = 0.7),
    makeIntegerParam(id = "min_split", lower=1L, upper=100L, default=15L),
    makeNumericParam(id = "min_percentage_split", lower=0, upper=1, default=0.25)
  ),
  properties = c(),
  run_fun = run_scuba,
  plot_fun = plot_scuba
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
  dir.create(temp_folder, recursive = TRUE)

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
        "cd {find.package('SCUBA')}/venv",
        "source bin/activate",
        "python3 {find.package('SCUBA')}/wrapper.py {temp_folder} 0 {c(0, 1)[as.numeric(rigorous_gap_stats)+1]} {N_dim} {low_gene_threshold} {low_gene_fraction_max} {min_split} {min_percentage_split}",
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
      to=as.character(to),
      directed=TRUE
    )

  milestone_ids <- unique(c(milestone_network$from, milestone_network$to))

  milestone_percentages <- data_frame(
    cell_id = rownames(counts),
    milestone_id = labels,
    percentage = 1
  )

  # remove temporary output
  unlink(temp_folder, recursive = TRUE)

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
  g <- ti_predictions$milestone_network %>% igraph::graph_from_data_frame()
  l <- igraph::layout_as_tree(g)
  igraph::plot.igraph(g, layout=l)
}
