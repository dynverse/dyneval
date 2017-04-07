library(SCORPIUS)
library(cowplot)
library(fitdistrplus)
library(tidyverse)
data(ginhoux)

expression <- ginhoux$expression
sample_info <- ginhoux$sample.info

dist <- correlation.distance(expression)
out_filt <- outlier.filter(dist)
dist <- dist[out_filt, out_filt]
expression <- expression[out_filt,]
sample_info <- sample_info[out_filt,,drop=F]

space <- reduce.dimensionality(dist)
traj <- infer.trajectory(space)

tapply(traj$time, sample_info$group.name, mean)
traj <- reverse.trajectory(traj)
time <- traj$time

draw.trajectory.plot(space, sample_info$group.name, traj$path)


state_names <- c("MDP", "CDP", "PreDC")

pct <- sapply(state_names, function(type) {
  typetime <- time[sample_info$group.name == type]
  fit <- fitdist(typetime,"norm")
  dnorm(time, fit$estimate[["mean"]], fit$estimate[["sd"]])
})

state_network <- tibble::data_frame(from = c("MDP", "CDP"), to = c("CDP", "PreDC"), length = 1)

max_froms <- lapply(state_names, function(sn_from) {
  sn_tos <- state_network %>% filter(from == sn_from) %>% .$to
  states <- unique(c(sn_from, sn_tos))
  list(states = states, values = rowSums(pct[,states,drop=F]))
})
max_from_mat <- sapply(max_froms, function(l) l$values)
which_max_froms <- apply(max_from_mat, 1, which.max)
ix <- bind_rows(lapply(names(which_max_froms), function(sample_name) {
  wm <- which_max_froms[[sample_name]]
  states <- max_froms[[wm]]$states
  notstates <- setdiff(state_names, states)
  data_frame(sample_name, notstates)
}))
pct[as.matrix(ix)] <- 0

pct <- pct / matrix(rep(rowSums(pct), each = 3), ncol = 3, byrow = T)

ggplot() + geom_line(aes(time, pct[,"MDP"], colour = "MDP")) + geom_line(aes(time, pct[,"CDP"], colour = "CDP")) + geom_line(aes(time, pct[,"PreDC"], colour = "PreDC"))


state_percentages <- data.frame(id = rownames(expression), pct) %>% as_data_frame

ginhoux_task <- wrap_ti_task_data(
  ti_type = "linear",
  name = "ginhoux",
  expression = expression,
  state_names = state_names,
  state_network = state_network,
  state_percentages = state_percentages,
  sample_info = sample_info
)
saveRDS(ginhoux_task, "data/ti/linear/ginhoux.rds")
