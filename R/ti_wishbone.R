#' Description for Wishbone
#' @export
description_wishbone <- function() create_description(
  name = "Wishbone",
  short_name = "Wishbone",
  package_loaded = c(),
  package_required = c("jsonlite", "Wishbone"),
  par_set = makeParamSet(
    makeIntegerParam(id = "knn", lower=2, upper=100, default=10),
    makeIntegerParam(id = "n_diffusion_components", lower=2, upper=20, default=10),
    makeIntegerParam(id = "n_pca_components", lower=2, upper=30, default=15),
    makeLogicalParam(id = "branch", default = TRUE),
    makeIntegerParam(id = "k", lower=2, upper=100, default=15),
    makeIntegerParam(id = "num_waypoints", lower=2, upper=500, default=250),
    makeLogicalParam(id = "normalize", default = TRUE),
    makeNumericParam(id = "epsilon", lower=0.1, upper=10, default=1)
  ),
  properties = c(),
  run_fun = run_wishbone,
  plot_fun = plot_wishbone
)

run_wishbone <- function(counts,
                         start_cell_id,
                         knn = 10,
                         n_diffusion_components = 2,
                         n_pca_components = 15,
                         markers="~",
                         branch=TRUE,
                         k=15,
                         num_waypoints=50,
                         normalize=TRUE,
                         epsilon=1
                      ) {
  if(is.null(start_cell_id)) stop("Give start cell id")

  requireNamespace("jsonlite")
  temp_folder <- tempfile()
  dir.create(temp_folder, recursive = T)

  counts %>%
    {log2(.+1)} %>%
    as.data.frame() %>%
    write.table(paste0(temp_folder, "/counts.tsv"), sep="\t")

  params <- as.list(environment())[formalArgs(run_wishbone)]
  params <- params[-(names(params) == "counts")]
  params[["components_list"]] <- seq_len(n_diffusion_components)-1
  params %>% jsonlite::toJSON(auto_unbox=TRUE) %>% write(paste0(temp_folder, "/params.json"))

  output <- system2(
    "/bin/bash",
    args = c(
      "-c",
      shQuote(glue::glue(
        "cd {find.package('Wishbone')}/venv",
        "source bin/activate",
        "python {find.package('Wishbone')}/wrapper.py {temp_folder}",
        .sep = ";"))
    ), stdout = T, stderr = T
  )

  print(output)

  branch_assignment <- jsonlite::read_json(paste0(temp_folder, "/branch.json")) %>%
    unlist() %>%
    {tibble(branch=., cell_id=names(.))}
  trajectory <- jsonlite::read_json(paste0(temp_folder, "/trajectory.json")) %>%
    unlist() %>%
    {tibble(time=., cell_id=names(.))}
  model <- left_join(branch_assignment, trajectory, by="cell_id")
  space <- read_csv(paste0(temp_folder, "/dm.csv")) %>% rename(cell_id=X1) %>% rename_if(is.numeric, funs(paste0("Comp", .)))

  if(branch) {
    milestone_network <- tibble(from=c("M1", "M2", "M2"), to=c("M2", "M3", "M4"), branch=c(1, 2, 3))
  } else {
    milestone_network <- tibble(from=c("M1"), to=c("M2"), branch=c(1))
  }

  progressions <- left_join(model, milestone_network, by="branch")

  # get lengths of milestone network
  milestone_network <- progressions %>% group_by(branch) %>% summarise(length=max(time) - min(time)) %>% left_join(milestone_network, by="branch")

  # now scale the times between 0 and 1 => percentages
  progressions <- progressions %>%
    group_by(branch) %>%
    mutate(percentage=(time - min(time))/(max(time) - min(time))) %>%
    ungroup()

  milestone_ids <- unique(c(milestone_network$from, milestone_network$to))
#
#   ggplot(left_join(model, space, by="cell_id")) + geom_point(aes(Comp0, Comp1, color=branch)) + viridis::scale_color_viridis()
#   ggplot(left_join(model, space, by="cell_id")) + geom_point(aes(Comp0, Comp1, color=time)) + viridis::scale_color_viridis()
#   ggplot(left_join(tasks$progressions[[1]], space, by="cell_id")) + geom_point(aes(Comp0, Comp1, color=percentage)) + viridis::scale_color_viridis()
#
#   plot(tasks$progressions[[1]]$percentage %>% unlist(), space$Comp1)

  # remove temporary output
  unlink(temp_folder, recursive = T)

  wrap_ti_prediction(
    ti_type = "linear",
    id = "Wishbone",
    cell_ids = rownames(counts),
    milestone_ids = milestone_ids,
    milestone_network = milestone_network %>% select(from, to, length),
    progressions = progressions %>% select(cell_id, from, to, percentage),
    space = space,
    model=model
  )
}

plot_wishbone <- function(ti_predictions) {
  ggplot(left_join(ti_predictions$model, ti_predictions$space, by="cell_id")) + geom_point(aes(Comp0, Comp1, color=branch))
}
