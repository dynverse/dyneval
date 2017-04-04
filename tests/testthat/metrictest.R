data_folder <- "tests/testthat/data_metrictest"
files <- list.files(data_folder)

trajectory_data <- lapply(files, function(name) {
  gold_path <- paste0(data_folder, "/", name, "/gold.ods")
  traj <- dyneval:::read_ods_trajectory(gold_path)
  dimensionality_reduce_trajectory(name, traj)
})

all_trajs <- list(
  milestones = bind_rows(trajectory_data %>% map(~ .$milestones)),
  lines = bind_rows(trajectory_data %>% map(~ .$lines)),
  cells = bind_rows(trajectory_data %>% map(~ .$cells))
)

make_trajectory_plot(trajectory[[8]])
with(trajectory_data[[8]], make_trajectory_plot(milestones, lines, cells))

# pdf("Rplot.pdf", 16, 12)
make_trajectory_plot(all_trajs) + facet_wrap(~name)
# dev.off()
