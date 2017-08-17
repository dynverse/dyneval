# summarise the scores
# get gs-toy score
# get difference and fractions between gs-toy and gs-gs
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
  geom_boxplot(aes(toy_category, diff, color=perturbator_id)) +
  facet_wrap(~score_id) +
  coord_flip()

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
  ggbeeswarm::geom_beeswarm(aes(toy_category, diff, color=perturbator_id)) +
  geom_hline(yintercept = 0) +
  facet_wrap(~score_id)

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

scores_summary %>%
  mutate(is_gs = (perturbator_id == "gs")) %>%
  ggplot() +
  geom_boxplot(aes(is_gs, diff, color=perturbator_id)) +
  facet_wrap(~score_id)

### Linear vs cycle
## 3a: Breaking of cycles should lower score
scores_summary %>%
  filter(perturbator_id == "break_cycles") %>%
  group_by(score_id) %>%
  summarise(rule_id="3a", rule = all(diff < 0)) %>%
  add_rule()

## 3b: Joing a linear trajectory to a cycle should lower score
scores_summary %>%
  filter(perturbator_id == "join_linear") %>%
  group_by(score_id) %>%
  summarise(rule_id="3b", rule = all(diff < 0)) %>%
  add_rule()

### 4: Large perturbations should have a larger decrease compared to a small perturbations
## 4a Shuffling cells
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
    rule_id = "4a",
    rule = all(maxdiff < 0)
  ) %>%
  add_rule()

scores_summary_largevssmall %>% ggplot() +
  geom_boxplot(aes(toy_category, score, color=perturbator_id)) +
  facet_wrap(~score_id)

## 4b Added states
scores_summary_largevssmall <- scores_summary %>%
  group_by(generator_id) %>%
  filter(perturbator_id %in% c("hairy_small", "hairy_large")) %>%
  filter(length(unique(perturbator_id)) == 2) %>%
  ungroup() %>%
  mutate(is_small = (perturbator_id == "hairy_small"))

scores_summary_largevssmall %>%
  group_by(score_id) %>%
  summarise(
    maxdiff = max(score[!is_small]) - min(score[is_small]),
    rule_id = "4b",
    rule = all(maxdiff < 0)
  ) %>%
  add_rule()

scores_summary_largevssmall %>% ggplot() +
  geom_boxplot(aes(toy_category, score, color=perturbator_id)) +
  facet_wrap(~score_id)

### 5: Splitting linear into bifurcating should lower score
scores_summary %>%
  filter(perturbator_id == "split_linear") %>%
  group_by(score_id) %>%
  summarise(rule_id="5", rule = all(diff < 0)) %>%
  add_rule()

### 6: Warping should lower score
scores_summary %>%
  filter(perturbator_id == "warp") %>%
  group_by(score_id) %>%
  summarise(rule_id="5", rule = all(diff < 0)) %>%
  add_rule()

scores_summary %>%
  filter(perturbator_id == "warp") %>%
  ggplot() +
  ggbeeswarm::geom_beeswarm(aes(score_id, diff, color=generator_id))


### 7: Changes in the milestone network should affect the scores, even if they do not greatly impact the ordering of the cells
## These are rules made to compare cell-based metrics with milestone network metrics

scores_summary %>%
  group_by(score_id) %>%
  mutate(maxscore = max(score)) %>%
  filter(perturbator_id == "hairy") %>%
  summarise(rule_id="7a", rule = all(score < maxscore/2)) %>%
  add_rule()






## Plot all rules
rules %>% ggplot(aes(rule_id, score_id)) +
  geom_raster(aes(fill = c("red", "green")[as.numeric(rule)+1])) +
  geom_text(aes(label = c("▼", "▲")[as.numeric(rule)+1], size=as.numeric(rule)), color="white") +
  scale_fill_identity()


## Possible score combinations satisfying all rules
allscores <-unique(rules$score_id)
score_ids_combinations <- map(1:length(allscores), function(x) map(combn(length(allscores), x, simplify = FALSE), ~allscores[.])) %>% unlist(recursive=FALSE)

map(score_ids_combinations, function(score_ids_combination) {
  rules %>% filter(score_id %in% score_ids_combination) %>%
    group_by(rule_id) %>%
    summarise(any=any(rule), number=sum(rule, na.rm=TRUE)) %>%
    group_by() %>%
    summarise(all=all(any), avg_rules_retrieved=sum(number)/length(score_ids_combination), nscores = length(score_ids_combination)) %>%
    mutate(score_ids_combination=list(score_ids_combination))
}) %>% bind_rows() %>% filter(all) %>% arrange(nscores, -avg_rules_retrieved) %>% View
