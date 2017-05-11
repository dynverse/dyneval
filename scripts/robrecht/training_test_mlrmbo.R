library(cowplot)
library(tidyverse)
library(dyneval)
library(mlrMBO)

datasets_info <- readRDS("../dyngen/results/datasets.rds")

output_root_folder <- "results/output_dyngentest/"
.datasets_location <- "../dyngen/results"

task_wrapped <- lapply(seq_len(nrow(datasets_info)), function(dataset_num) {
  dataset_id <- datasets_info$id[[dataset_num]]
  dataset <- dyngen::load_dataset(dataset_id, contents = dyngen::contents_dataset(experiment = F))

  with(dataset, dyneval::wrap_ti_task_data(
    ti_type = model$modulenetname,
    name = info$id,
    counts = counts,
    state_names = gs$milestone_names,
    state_net = gs$milestone_net,
    state_percentages = gs$milestone_percentages %>% slice(match(rownames(counts), id))
  ))
})
make_obj_fun <- function(method) {
  makeSingleObjectiveFunction(
    name = "TItrain",
    vectorized = F,
    minimize = F,
    noisy = T,
    has.simple.signature = F,
    par.set = method$par_set,
    fn = function(x, tasks) {
      outs <- lapply(tasks, function(task) {
        arglist <- c(list(counts = task$counts), x)
        model <- do.call(method$run_fun, arglist)
        coranking <- compute_coranking(task$geodesic_dist, model$geodesic_dist)
        list(model = model, coranking = coranking)
      })
      outs %>% map_dbl(~ .$coranking$summary$max_lcmc) %>% mean
    })
}

method <- description_scorpius()
control <- makeMBOControl(propose.points = 8)
design <- generateDesign(n = 5*8, par.set = method$par_set)
tasks <- task_wrapped[c(1,2,8)]

parallelStartMulticore(cpus = 8, show.info = TRUE)
tune_out <- mbo(make_obj_fun(method), design = design, control = control, show.info = T, more.args = list(tasks = tasks))
parallelStop()

tune_out



