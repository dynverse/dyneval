library(cowplot)
library(mlr)
library(tidyverse)
library(dyneval)

datasets_info <- readRDS("../dyngen/results/datasets.rds")

output_root_folder <- "results/output_dyngentest/"
.datasets_location <- "../dyngen/results"

tasks <- lapply(seq_len(nrow(datasets_info)), function(dataset_num) {
  dataset_id <- datasets_info$id[[dataset_num]]
  dataset <- dyngen::load_dataset(dataset_id, contents = dyngen::contents_dataset(experiment = F))

  expression <- as.data.frame(log2(dataset$counts+1))
  gold <- list(
    state_names = dataset$gs$milestone_names,
    state_network = dataset$gs$milestone_net,
    state_percentages = dataset$gs$milestone_percentages %>% slice(match(rownames(dataset$counts), id))
  )

  makeTITask(id = paste0(dataset$model$modulenetname, "_", dataset$info$id), data = expression, gold = gold)
})

task_wrapped <- lapply(seq_len(nrow(datasets_info)), function(dataset_num) {
  dataset_id <- datasets_info$id[[dataset_num]]
  dataset <- dyngen::load_dataset(dataset_id, contents = dyngen::contents_dataset(experiment = F))

  with(dataset, dyneval::wrap_ti_task_data(
    ti_type = model$modulenetname,
    name = info$id,
    expression = log2(counts+1),
    state_names = gs$milestone_names,
    state_net = gs$milestone_net,
    state_percentages = gs$milestone_percentages %>% slice(match(rownames(counts), id))
  ))
})

# try with scorpius and the first task
learner <- makeLearner("ti.scorpius")
task <- tasks[[1]]
task_wrap <- task_wrapped[[1]]


# train using default params
model <- train(learner, task)
performance(pred = NULL, measures = list(corank_auc, corank_max, geodesic_cor), task = task, model = model)

# train the parameters in a grid
ctrl <- makeTuneControlGrid(resolution = 5L)
rdesc <- makeResampleDesc("CV", iters = 5L)
res <- tuneParams(learner, task, rdesc, par.set = learner$par.set, control = ctrl)

# train a new model using these params
par_values <- res$x
model <- train(makeLearner("ti.scorpius", par.vals = par_values), task)

g0 <- plotLearner.ti.default(task_wrap)
g1 <- plotLearner.ti.default(model$learner.model)
g2 <- plotLearner.ti.scorpius(model$learner.model)
g3 <- plotLearner.ti.combined(task_wrap, model$learner.model)
cowplot::plot_grid(g0, g1, g2, g3)




# benchmarking?
learner <- makeLearner("ti.scorpius")
rdesc <- makeResampleDesc("CV", iters = 2L)
bm <- benchmark(learner, tasks, rdesc)


ctrl = makeTuneControlIrace(maxExperiments = 200)
ctrl = makeTuneControlGrid(resolution = 10L)
# ctrl = makeTuneControlMBO()
rdesc = makeResampleDesc("CV", iters = 2L)
# rdesc = makeResampleDesc("Bootstrap", iters = 2L)
res = tuneParams(learner, task, rdesc, par.set = learner$par.set, control = ctrl)
