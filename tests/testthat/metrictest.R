library(tidyverse)
library(dyneval)

data_folder <- "tests/testthat/data_metrictest"
files <- list.files(data_folder)

trajectory_data <- lapply(files, function(name) {
  cat("Reading ", name, "\n", sep = "")
  gold_path <- paste0(data_folder, "/", name, "/gold.ods")
  traj <- dyneval:::read_ods_trajectory(gold_path)
  dimensionality_reduce_trajectory(name, traj)
})

all_trajs <- list(
  milestones = bind_rows(trajectory_data %>% map(~ .$milestones)),
  lines = bind_rows(trajectory_data %>% map(~ .$lines)),
  cells = bind_rows(trajectory_data %>% map(~ .$cells))
)

# make_trajectory_plot(trajectory_data[[7]])
# make_trajectory_plot(trajectory_data[[6]])

# pdf("alltestplots.pdf", 16, 12)
make_trajectory_plot(all_trajs) + facet_wrap(~name, ncol = 5)
# dev.off()


gold_path <- paste0(data_folder, "/1a_linear/gold.ods")
pred_path <- paste0(data_folder, "/1a_linear/prediction_good.ods")
traj <- dyneval:::read_ods_trajectory(gold_path)
pred_traj <- dyneval:::read_ods_trajectory(pred_path)
dimensionality_reduce_trajectory(name, traj)

# # test renaming
# traj <- list(
#   structure = data.frame(from = "x", to = "y", length = 1, stringsAsFactors = F),
#   cells = data_frame(x = seq(0, 1, by = .1), y = 1 - x)
# )
# drtraj <- dimensionality_reduce_trajectory("test", traj)
# make_trajectory_plot(drtraj)
#
# # test more than 12 milestones
# traj <- list(
#   structure = data.frame(from = letters[-26], to = letters[-1], length = 1, stringsAsFactors = F),
#   cells = local({
#     d <- diag(length(letters))
#     colnames(d) <- letters
#     as.data.frame(d)
#   })
# )
# drtraj <- dimensionality_reduce_trajectory("test", traj)
# make_trajectory_plot(drtraj)
