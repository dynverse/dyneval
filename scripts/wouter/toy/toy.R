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
  "bifurcating", "gs",
  "bifurcating", "switch_two_cells",
  "bifurcating", "switch_all_cells",
  "cycle", "gs",
  "cycle", "switch_all_cells",
  "cycle", "break_cycles"
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

# generate gold standards and toys, can take some time (for computing the geodesic distances)
# I choose to not do this using mutate because it is much easier to debug
toys$gs <- map2(toys$generator, toys$ncells, ~.x(.y))
toys$toy <- map2(toys$perturbator, toys$gs, ~.x(.y))
toys$toy <- map2(toys$toy, toys$toy_id, ~rename_toy(.x, .y))

# plot toys
toys <- toys %>% rowwise() %>% mutate(plot_strip=list(plot_strip(gs, toy)))
cowplot::plot_grid(plotlist=toys$plot_strip, ncol=nreplicates)

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
  fun <- make_obj_fun(dummy_method(gs))
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

# summarise the scores
# get gs-toy score
# get difference and fractions between gs-toy and gs-gs

scores <- toys %>% rowwise() %>% do(compare_gs_toy(.$gs, .$toy))
scores_summary <- scores %>% gather(score_id, score, -toy_id, -comparison) %>%
  spread(comparison, score) %>% mutate(
  score=`gs-toy`,
  diff= `gs-toy`-`gs-gs`,
  frac=`gs-toy`/`gs-gs`
) %>% select(-`gs-toy`, -`gs-gs`) %>% left_join(toys, by="toy_id")

scores_summary %>%
  ggplot() +
    geom_boxplot(aes(toy_category, score, color=perturbator_id)) +
    facet_wrap(~score_id) +
    coord_flip()

scores_summary %>%
  ggplot() +
  geom_boxplot(aes(toy_category, diff)) +
  facet_wrap(~score_id)

rules <- tibble()
add_rule <- function(rule) {
  rules <<- bind_rows(rules, select(rule, score_id, rule_id, rule)) %>%
    group_by(rule_id, score_id) %>%
    filter(row_number() == n()) %>%
    ungroup()
  rule
}

### A good score...
### 1: Same score on any gold standard, irrespective of structure or number of cells
equal_tol <- function(a, b, tol = 1e-5) {abs(a-b) < tol}

scores_summary %>%
  filter(perturbator_id == "gs") %>%
  group_by(score_id) %>%
  summarise(rule_id = "1", rule = equal_tol(max(score), min(score))) %>%
  add_rule()

scores_summary %>%
  filter(perturbator_id == "gs") %>%
  ggplot() +
  ggbeeswarm::geom_beeswarm(aes(score_id, score, color=toy_category))

### Gold standard versus prediction
## 2a: Score should be lower before and after perturbation (direct gs comparison)
scores_summary %>%
  filter(perturbator_id != "gs") %>%
  group_by(score_id) %>%
  summarise(rule_id = "2a", rule = all(diff < 0)) %>%
  add_rule()

scores_summary %>%
  ggplot() +
    ggbeeswarm::geom_beeswarm(aes(score_id, diff, color=toy_category))

## 2b: Score should be lower before and after perturbation, irrespective of structure or number of cells (indirect gs comparison)
scores_summary %>%
  mutate(is_gs = (perturbator_id == "gs")) %>%
  group_by(score_id) %>%
  summarise(
    maxdiff = max(score[!is_gs]) - min(score[is_gs]),
    rule_id = "2b",
    rule = maxdiff <= 0
  ) %>%
  add_rule()

scores_summary %>%
  mutate(is_gs = (perturbator_id == "gs")) %>%
  ggplot() +
  geom_boxplot(aes(is_gs, diff, color=toy_category)) +
  facet_wrap(~score_id)

## Linear vs cycle
## 3a: Breaking of cycles should lower score
scores_summary %>%
  filter(perturbator_id == "break_cycles") %>%
  group_by(score_id) %>%
  summarise(rule_id="3", rule = all(diff < 0)) %>%
  add_rule()

## 3b: Joing a linear trajectory to a cycle should lower score


## 4: Large perturbations should have a larger decrease compared to a small perturbations
scores_summary_largevssmall <- scores_summary %>%
  group_by(generator_id) %>%
  filter(perturbator_id %in% c("switch_two_cells", "switch_all_cells")) %>%
  filter(length(unique(perturbator_id)) == 2) %>%
  ungroup() %>%
  mutate(is_small = (perturbator_id == "switch_two_cells"))

scores_summary_largevssmall %>%
  group_by(score_id) %>%
  summarise(
    maxdiff = max(score[!is_small]) - min(score[is_small]),
    rule_id = "4",
    rule = all(maxdiff < 0)
  ) %>%
  add_rule()

scores_summary_largevssmall %>% ggplot() +
  geom_boxplot(aes(toy_category, score, color=perturbator_id)) +
  facet_wrap(~score_id)

##
rules %>% ggplot(aes(rule_id, score_id)) +
  geom_raster(aes(fill = c("red", "green")[as.numeric(rule)+1])) +
  geom_text(aes(label = c("✘", "✔")[as.numeric(rule)+1]), color="white") +
  scale_fill_identity()
