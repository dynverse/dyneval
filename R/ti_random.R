#' Description for shuffled
#' @export
description_shuffled <- function() create_description(
  name = "shuffle",
  short_name = "shuffle",
  package_loaded = c(),
  package_required = c(),
  par_set = makeParamSet(
    makeIntegerParam(id = "dummy_param", lower = 0, default = 0.5, upper = 1)
  ),
  properties = c(),
  run_fun = run_shuffled,
  plot_fun = plot_default
)

run_shuffled <- function(counts, task, dummy_param = .5) {
  allcells <- rownames(counts)
  mapper <- setNames(allcells, sample(allcells))
  progressions <- task$progressions
  progressions$cell_id <- mapper[progressions$cell_id]

  wrap_ti_prediction(
    ti_type = task$ti_type,
    id = paste0(task$id, "_shuffled"),
    cell_ids = task$cell_ids,
    milestone_ids =task$milestone_ids,
    milestone_network = task$milestone_network,
    progressions = progressions
  )
}

