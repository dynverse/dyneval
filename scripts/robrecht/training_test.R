library(cowplot)
library(mlr)
library(tidyverse)
library(dyneval)

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
task_tib <- to_tibble(task_wrapped)

task <- makeTITask(id = "dyngen", data = task_tib[c(1,6,9,20),])

learner <- makeLearner("ti.scorpius")
# learner <- makeLearner("ti.monocleDDRtree")
# learner <- makeLearner("ti.tscan")

# # train using default params
# model <- train(learner, task)
# pred <- predict(model, task)
#
# performance(pred, measures = list(corank_auc, corank_max, geodesic_cor))

# train the parameters in a grid
# ctrl <- makeTuneControlGrid(resolution = 2L)
ctrl <- makeTuneControlMBO()
rdesc <- makeResampleDesc("CV", iters = 2L)
par <- makeParamSet(
  makeIntegerParam(id = "num_dimensions", lower = 2L, default = 3L, upper = 20L),
  makeIntegerParam(id = "num_clusters", lower = 2L, default = 4L, upper = 20L),
  makeDiscreteParam(id = "distance_method", default = "spearman", values = c("spearman", "pearson", "kendall"))
)

configureMlr(on.learner.warning = "quiet", show.learner.output = FALSE)

library(parallelMap)
# parallelStartMulticore(cpus = 8, show.info = TRUE)
res <- tuneParams(learner = learner, task = task, resampling = rdesc, measures = corank_max, par.set = par, control = ctrl)
# parallelStop()









