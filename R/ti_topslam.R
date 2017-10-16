#' Description for topslam
#' @export
description_topslam <- function() create_description(
  name = "topslam",
  short_name = "topslam",
  package_loaded = c(),
  package_required = c("jsonlite", "topslam"),
  par_set = makeParamSet(
    makeIntegerParam(id = "n_components", lower=2L, upper=10L, default=2L),
    makeIntegerParam(id = "n_neighbors", lower=2L, upper=100L, default=10L),
    makeIntegerParam(id = "linear_dims", lower=0L, upper=5L, default=0L),
    makeIntegerParam(id = "max_iters", lower=100L, upper=1000L, default=200L),
    makeLogicalVectorParam(id = "dimreds", len = 5, default = rep(TRUE, 5))
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
                        max_iters = 200,
                        dimreds = rep(TRUE, 5)
                      ) {
  requireNamespace("jsonlite")
  temp_folder <- tempfile()
  dir.create(temp_folder, recursive = TRUE)

  counts %>%
    {log2(.+1)} %>%
    as.data.frame() %>%
    write.table(paste0(temp_folder, "/counts.tsv"), sep="\t")

  params <- as.list(environment())[formalArgs(run_topslam)]
  params <- params[-(names(params) == "counts")]
  params$start_cell_id <- which(rownames(counts) == start_cell_id) -1 # add this minus 1 because R starts counting from 1
  params$dimreds <- c("t-SNE", "PCA", "Spectral", "Isomap", "ICA")[dimreds]
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
    ), stdout = TRUE, stderr = TRUE
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
  unlink(temp_folder, recursive = TRUE)

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

#' @importFrom viridis scale_colour_viridis
plot_topslam <- function(prediction) {
  ggplot() +
    geom_raster(aes(x, y, fill=energy), data=prediction$wad) +
    geom_contour(aes(x, y, z=energy, weight=energy), data=prediction$wad, binwidth = 0.05, color="#222222", alpha=0.4) +
    geom_point(aes(Comp1, Comp2, color=time), data=prediction$model) +
    scale_fill_gradientn(colors=c("grey75", "white")) +
    viridis::scale_colour_viridis()
}
