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
    package_installed = c("jsonlite"),
    par_set = makeParamSet(
## TODO
    ),
    properties = c(),
    run_fun = run_scuba,
    plot_fun = plot_scuba
  )
}

#' @importFrom utils write.table
run_scuba <- function(counts,
                      log_mode=TRUE,
                      rigorous_gap_stats=TRUE) {

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
        "python {path.package('dyneval')}/extra_code/PySCUBA/scuba_wrapper.py {temp_folder} {log_mode} {rigorous_gap_stats}",
        .sep = ";"))
    )
  )

  output <- jsonlite::read_json(paste0(temp_folder, "/output.json"))
  labels <- as.character(unlist(output$labels))
  tree <- output$tree

  milestone_network <- map2(names(tree), tree, ~tibble(from=.x, to=.y)) %>%
    bind_rows() %>%
    unnest(to) %>%
    filter(to > -1) %>%
    mutate(to=as.character(to)) %>%
    mutate(length=1)

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
