source("../dyngen/scripts/evaluation/methods/dimred.R")
library(tidyverse)

expr <- read_tsv("/home/wouters/thesis/projects/dyneval/data/real/fibroblast_reprogramming/GSE67310_iN_data_log2FPKM_annotated.txt")
expression <- expr[, -c(1:5)] %>% as.matrix() %>% magrittr::set_rownames(expr$cell_name)
cellinfo <- expr[, c(1:5)] %>% as.data.frame() %>% magrittr::set_rownames(expr$cell_name)

datasetinfo <- list(organism = "mouse", genenames = "symbol", technology = "fluidigm_c1", name="fibroblast_reprogramming_treutlein")

dataset <- lst(expression, cellinfo, datasetinfo)

cellinfo$assignment %>% table


statenet = tribble(
  ~from, ~to,
  "MEF", "d2_intermediate",
  "d2_intermediate", "d2_induced",
  "d2_induced", "d5_intermediate",
  "d5_intermediate", "d5_earlyiN",
  "d5_earlyiN", "Neuron",
  "d5_earlyiN", "Myocyte"
)





cellsoi = cellinfo %>% filter((assignment %in% c("MEF", "d2_induced", "d2_intermediate", "d5_earlyiN", "d5_intermediate", "Myocyte", "Neuron"))) %>% pull(cell_name)
space = SCORPIUS::reduce.dimensionality(SCORPIUS::correlation_distance(expression[cellsoi, ]), ndim=3)
space = tsne(expression[cellsoi, ], ndim=2)
traj = SCORPIUS::infer.trajectory(space)
#space = ica(expression)
SCORPIUS::draw.trajectory.plot(space, cellinfo[cellsoi, ]$assignment, traj$path)

assignment_int <- as.numeric(factor(cellinfo[cellsoi, ]$assignment))

library(rgl)
plot3d(space[, 1], space[, 2], space[, 3], col=rainbow(max(assignment_int))[assignment_int])
