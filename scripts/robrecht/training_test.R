library(mlr)
library(tidyverse)
library(dyneval)

datasets_info <- readRDS("../dyngen/results/datasets.rds")

output_root_folder <- "results/output_dyngentest/"
.datasets_location <- "../dyngen/results"

tasks <- lapply(seq_len(nrow(datasets_info)), function(dataset_num) {
  dataset_id <- datasets_info$id[[dataset_num]]
  dataset <- dyngen::load_dataset(dataset_id, contents = dyngen::contents_dataset(experiment = F))

  task <- with(dataset, dyneval::wrap_ti_task_data(
    ti_type = model$modulenetname,
    name = info$id,
    expression = log2(counts+1),
    state_names = gs$milestone_names,
    state_net = gs$milestone_net,
    state_percentages = gs$milestone_percentages %>% slice(match(rownames(counts), id))
  ))
})


learner <- makeLearner("ti.scorpius")
task <- makeTITask(data = as.data.frame(tasks[[1]]$expression))
mod <- train(learner, task)

print(mod)
## use random forest to classify iris data
task = makeClassifTask(data = iris, target = "Species")
learner = makeLearner("classif.rpart", minsplit = 7, predict.type = "prob")
mod = train(learner, task, subset = training.set)
