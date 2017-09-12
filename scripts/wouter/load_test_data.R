library(tidyverse)
library(dyneval)

check_dependencies()

set.seed(1)
tasks <- generate_toy_datasets()

i <- 2
counts <- tasks$counts[[i]]
special_cells <- tasks$special_cells[[i]]
cell_grouping <- tasks$cell_grouping[[i]]
progressions <- tasks$progressions[[i]]
