#!/usr/local/bin/Rscript

library(optparse)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(dyneval)

metrics <- dyneval::metrics

parser <-
  OptionParser(usage = paste0("LOCAL=/path/to/folder; MOUNT=/ti; docker run -v $LOCAL:$MOUNT dynverse/dyneval")) %>%
  add_option("--dataset", type = "character", help = "Filename of the dataset, example: $MOUNT/dataset.(h5|loom). h5 files can be created by cellranger or dyncli.") %>%
  add_option("--model", type = "character", help = "Filename of the model, example: $MOUNT/dataset.(h5|loom). h5 files can be created by cellranger or dyncli.")  %>%
  add_option("--output", type = "character", help = "Filename of the scores, example: $MOUNT/dataset.(h5|loom). Will be a json file containing the scores.") %>%
  add_option("--metrics", type = "character", help = "Which metrics to calculate, example: correlation,him,F1_milestones,featureimp_wcor", default = "correlation,him,F1_milestones,featureimp_wcor")

parsed_args <- parse_args(parser, args = commandArgs(trailingOnly = TRUE))

if (any(sapply(parsed_args[c("dataset", "model", "output")], is.null))) {
  stop("dataset, model and output arguments are mandatory")
}

# read dataset and model
dataset <- dynutils::read_h5(parsed_args$dataset)
model <- dynutils::read_h5(parsed_args$model)

metrics <- strsplit(parsed_args$metrics, ",")[[1]]

scores <- dyneval::calculate_metrics(dataset, model, metrics = metrics)

# write output
jsonlite::write_json(scores, parsed_args$output)
