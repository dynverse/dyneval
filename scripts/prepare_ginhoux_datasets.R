library(SCORPIUS)
library(cowplot)
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

traj <- reverse.trajectory(traj)
time <- traj$time

library(fitdistrplus)
pct <- sapply(c("MDP", "CDP", "PreDC"), function(type) {
  typetime <- time[sample_info$group.name == type]
  fit <- fitdist(typetime,"norm")
  dnorm(time, fit$estimate[["mean"]], fit$estimate[["sd"]])
})
pct <- pct / matrix(rep(rowSums(pct), each = 3), ncol = 3, byrow = T)

ggplot() + geom_line(aes(time, pct[,"MDP"], colour = "MDP")) + geom_line(aes(time, pct[,"CDP"], colour = "CDP")) + geom_line(aes(time, pct[,"PreDC"], colour = "PreDC"))

state_network <- tibble::data_frame(from = c("MDP", "CDP"), to = c("CDP", "PreDC"), length = 1)
state_percentages <- data.frame(cellid = rownames(expression), pct) %>% as_data_frame

ginhoux_task <- dyneval:::make_ti_task_data("linear", "ginhoux", expression, state_network, state_percentages, sample_info)
saveRDS(ginhoux_task, "data/linear_ginhoux.rds")
