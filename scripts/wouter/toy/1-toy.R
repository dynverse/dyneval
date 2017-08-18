library(tidyverse)
library(dyneval)

source("scripts/wouter/toy/generation.R")
source("scripts/wouter/toy/perturbation.R")
source("scripts/wouter/toy/plotting.R")

toys_blueprint <- tribble(
  ~generator_id, ~perturbator_id,
  "linear", "gs",
  "linear", "switch_two_cells",
  "linear", "switch_all_cells",
  "linear", "join_linear",
  "linear", "split_linear",
  "linear", "warp",
  "linear", "hairy",
  "linear", "hairy_small",
  "linear", "hairy_large",
  "bifurcating", "gs",
  "bifurcating", "switch_two_cells",
  "bifurcating", "switch_all_cells",
  "bifurcating", "warp",
  "bifurcating", "hairy",
  "cycle", "gs",
  "cycle", "switch_all_cells",
  "cycle", "break_cycles",
  "cycle", "warp",
  "cycle", "hairy"
) %>% rowwise() %>% mutate(
  generator=list(get(paste0("generate_", generator_id))),
  perturbator=list(get(paste0("perturb_", perturbator_id)))
)

# replicate
nreplicates <- 5
toys <- toys_blueprint %>% slice(rep(1:n(), each=nreplicates)) %>% mutate(
  replicate=seq_len(nrow(.))%%nreplicates,
  toy_category=paste0(generator_id, "-", perturbator_id),
  toy_id=paste0(toy_category, "-", replicate)
) %>% mutate(ncells=runif(n(), 10, 200))

# generate gold standards and toys, can take some time (for computing the geodesic distances I presume)
# I choose to not do this using mutate because it is much easier to debug
toys$gs <- toys %>% split(seq_len(nrow(toys))) %>% parallel::mclapply(function(row) {
  row$generator[[1]](row$ncells)
}, mc.cores=8)
#toys$gs <- map2(toys$generator, toys$ncells, ~.x(.y))

toys$toy <- toys %>% split(seq_len(nrow(toys))) %>% parallel::mclapply(function(row) {
  row$perturbator[[1]](row$gs[[1]])
}, mc.cores=8)
#toys$toy <- map2(toys$perturbator, toys$gs, ~.x(.y))
toys$toy <- map2(toys$toy, toys$toy_id, ~rename_toy(.x, .y))

# plot toys
toys <- toys %>% rowwise() %>% mutate(plot_strip=list(plot_strip(gs, toy))) %>% ungroup()
cowplot::plot_grid(
  plotlist=toys %>% group_by(toy_category) %>% filter(row_number()==1) %>% pull(plot_strip)
)

# get the scores when comparing the gs to toy
compare_toy <- function(gs, toy, id=toy$id) {
  dummy_method <- function(pred) list(
    name = "dummy",
    short_name = "dummy",
    package_load = c(),
    package_installed = c(),
    par_set = ParamHelpers::makeParamSet(),
    properties = c(),
    run_fun = function(counts) pred,
    plot_fun = function(task) plot(1:10)
  )
  toy$counts <- "dummy"
  fun <- make_obj_fun(dummy_method(gs), metrics = c("Q_global", "Q_local", "correlation"))
  result <- fun(list(), dyneval:::list_as_tibble(list(toy)))
  scores <- attr(result, "extras")$.summary
  scores %>% select(-task_id) %>% mutate(toy_id=id)
}

# get the scores when both comparing gs to gs and comparing gs to toy
compare_gs_toy <- function(gs, toy) {
  bind_rows(
    compare_toy(gs, toy) %>% mutate(comparison="gs-toy"),
    compare_toy(gs, gs, id=toy$id) %>% mutate(comparison="gs-gs")
  )
}

scores <- toys %>% rowwise() %>% do(compare_gs_toy(.$gs, .$toy))














scores_summary %>% filter(is.na(score)) %>% pull(toy_id)
