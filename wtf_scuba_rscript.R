library(tidyverse)
library(dyneval)

tasks <- generate_toy_datasets()
tasks <- tasks[1,]

counts <- tasks$counts[[1]]
special_cells <- tasks$special_cells[[1]]
progressions <- tasks$progressions[[1]]

# choose certain parameters for each method, at which we know this method will perform well for the toy dataset
method_descriptions <- list(
  waterfall=list(),
  scorpius=list(),
  slingshot=list(),
  gpfates=list(nfates=1),
  stemid=list(clustnr=10, bootnr=10, pdishuf=10),
  tscan=list(),
  embeddr=list(nn_pct = 2),
  celltree_gibbs=list(sd_filter = 0),
  celltree_maptpx=list(sd_filter = 0),
  celltree_vem=list(sd_filter = 0),
  scuba=list(),
  slicer = list(min_branch_len=50),
  monocle_ddrtree=list(),
  wishbone = list(branch=F),
  random_linear=list()
)

method_descriptions <- method_descriptions["scuba"]

metric_names <- c("mean_R_nx", "auc_R_nx", "Q_local", "Q_global", "correlation", "isomorphic", "robbie_network_score")


dyneval:::run_scuba(counts)
