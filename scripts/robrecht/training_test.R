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

task <- makeTITask(id = "dyngen", data = task_tib)

learner <- makeLearner("ti.scorpius")
learner <- makeLearner("ti.monocleDDRtree")
learner <- makeLearner("ti.tscan")

# train using default params
model <- train(learner, task)
pred <- predict(model, task)

performance(pred, measures = list(corank_auc, corank_max, geodesic_cor))

# train the parameters in a grid
ctrl <- makeTuneControlGrid(resolution = 2L)
rdesc <- makeResampleDesc("CV", iters = 2L)
res <- tuneParams(learner, task, rdesc, par.set = learner$par.set, control = ctrl)

# # train a new model using these params
# par_values <- res$x
# model <- train(makeLearner("ti.scorpius", par.vals = par_values), task)
#
# g0 <- plotLearner.ti.default(task_wrap)
# g1 <- plotLearner.ti.default(model$learner.model)
# g2 <- plotLearner.ti.scorpius(model$learner.model)
# g3 <- plotLearner.ti.combined(task_wrap, model$learner.model)
# cowplot::plot_grid(g0, g1, g2, g3)
#
#
#
#
# ctrl = makeTuneControlIrace(maxExperiments = 200)
# ctrl = makeTuneControlGrid(resolution = 10L)
# # ctrl = makeTuneControlMBO()
# rdesc = makeResampleDesc("CV", iters = 2L)
# # rdesc = makeResampleDesc("Bootstrap", iters = 2L)
# res = tuneParams(learner, task, rdesc, par.set = learner$par.set, control = ctrl)
#
#
