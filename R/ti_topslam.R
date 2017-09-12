#' Description for topslam
#' @export
description_topslam <- function() create_description(
  name = "topslam",
  short_name = "topslam",
  package_loaded = c(),
  package_required = c("jsonlite", "topslam"),
  par_set = makeParamSet(
    # makeIntegerParam(id = "knn", lower=2, upper=100, default=10),
    # makeIntegerParam(id = "n_diffusion_components", lower=2, upper=20, default=10),
    # makeIntegerParam(id = "n_pca_components", lower=2, upper=30, default=15),
    # makeLogicalParam(id = "branch", default = TRUE),
    # makeIntegerParam(id = "k", lower=2, upper=100, default=15),
    # makeIntegerParam(id = "num_waypoints", lower=2, upper=500, default=250),
    # makeLogicalParam(id = "normalize", default = TRUE),
    # makeNumericParam(id = "epsilon", lower=0.1, upper=10, default=1)
  ),
  properties = c(),
  run_fun = run_topslam,
  plot_fun = plot_topslam
)

run_topslam <- function(counts,
                        start_cell_id,
                        n_components = 2,
                        n_neighbors = 10,
                        linear_dims = 0,
                        max_iters = 200
                      ) {
  requireNamespace("jsonlite")
  temp_folder <- tempfile()
  dir.create(temp_folder, recursive = T)

  counts %>%
    {log2(.+1)} %>%
    as.data.frame() %>%
    write.table(paste0(temp_folder, "/counts.tsv"), sep="\t")

  params <- as.list(environment())[formalArgs(run_topslam)]
  params <- params[-(names(params) == "counts")]
  params$start_cell_id <- which(rownames(counts) == start_cell_id) -1 # add this minus 1 because R starts counting from 1
  params %>% jsonlite::toJSON(auto_unbox=TRUE) %>% write(paste0(temp_folder, "/params.json"))

  output <- system2(
    "/bin/bash",
    args = c(
      "-c",
      shQuote(glue::glue(
        "cd {find.package('topslam')}/venv",
        "source bin/activate",
        "python {find.package('topslam')}/wrapper.py {temp_folder}",
        .sep = ";"))
    ), stdout = T, stderr = T
  )

  print(output)

  wad_grid <- read_csv(paste0(temp_folder, "/wad_grid.csv"))
  wad_energy <- read_csv(paste0(temp_folder, "/wad_energy.csv"))
  space <- read_csv(paste0(temp_folder, "/space.csv"))
  pseudotime <- read_csv(paste0(temp_folder, "/pseudotime.csv"))

  wad <- bind_cols(wad_grid, wad_energy)
  model <- bind_cols(space, pseudotime, tibble(cell_id = rownames(counts)))

  milestone_network <- tibble(from=c("M1"), to=c("M2"), length=1, directed=TRUE)
  progressions <- tibble(from="M1", to="M2", percentage=pseudotime$time, cell_id=rownames(counts))

  milestone_ids <- unique(c(milestone_network$from, milestone_network$to))

  # remove temporary output
  unlink(temp_folder, recursive = T)

  prediction <- wrap_ti_prediction(
    ti_type = "linear",
    id = "topslam",
    cell_ids = rownames(counts),
    milestone_ids = milestone_ids,
    milestone_network = milestone_network %>% select(from, to, length, directed),
    progressions = progressions %>% select(cell_id, from, to, percentage),
    model=model,
    wad=wad
  )
  prediction
}

plot_topslam <- function(prediction) {
  ggplot() +
    geom_raster(aes(x, y, fill=energy), data=prediction$wad) +
    geom_contour(aes(x, y, z=energy, weight=energy), data=prediction$wad, binwidth = 0.05, color="#222222", alpha=0.4) +
    geom_point(aes(Comp1, Comp2, color=time), data=model) +
    scale_fill_gradientn(colors=c("grey75", "white")) +
    viridis::scale_color_viridis()
}
