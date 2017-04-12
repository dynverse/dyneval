library(tidyverse)
library(dyneval)

my_read_ods <- function(path = NULL, sheet = 1, col_names = TRUE, col_types = NULL, na = "", skip = 0, formula_as_formula = FALSE, range = NULL) {
  res <- readODS:::parse_ods_to_sheets(path)
  ods_ns <- res[[2]]
  sheets <- res[[1]]
  target_sheet <- readODS:::select_sheet(sheets, ods_ns = ods_ns, which_sheet = sheet)
  cell_values <- readODS:::parse_rows(target_sheet, ods_ns, formula_as_formula = formula_as_formula, skip = skip)
  parsed_df <- readODS:::to_data_frame(cell_values = cell_values, header = col_names, na = na)
  if (length(col_types) == 1 && is.null(col_types)) {
    raw_sheet <- readr::type_convert(df = parsed_df)
  } else if (length(col_types) == 1 && is.na(col_types)) {
    raw_sheet <- parsed_df
  } else {
    raw_sheet <- readr::type_convert(df = parsed_df, col_types = col_types)
  }
  if (!is.null(range)) {
    res <- readODS:::select_range(raw_sheet, range)
  } else {
    res <- raw_sheet
  }
  return(res)
}

read_trajectory_ods <- function(name) {
  cat("Reading ", name, "\n", sep = "")
  filename <- paste0(data_folder, "/", name, "/gold.ods")

  structure_coltypes <- cols_only(from = col_character(), to = col_character(), length = col_double())
  cells_coltypes <- cols(cellid = col_character(), .default = col_double())

  # read structure tab
  structure <- my_read_ods(
    filename,
    sheet = "structure",
    col_types = structure_coltypes
  ) %>% mutate(
    from = ifelse(grepl("^[0-9]*$", from), paste0("milestone_", from), from),
    to = ifelse(grepl("^[0-9]*$", to), paste0("milestone_", to), to)
  )

  milestones <- sort(unique(c(structure$from, structure$to)))

  # read cells tab
  cells <- my_read_ods(
    filename,
    sheet = "cells",
    col_types = cells_coltypes
  ) %>%
    mutate(id = paste0("cell_", cellid)) %>%
    select(id, one_of(milestones))

  # fill in NAs
  for (milestone in milestones) {
    cells[,milestone] <- ifelse(is.na(cells[,milestone]), 0, cells[,milestone])
  }

  dyneval::wrap_ti_prediction(
    ti_type = "*",
    name = name,
    state_names = milestones,
    state_network = structure,
    state_percentages = cells,
    task_id = NULL
  )
}

data_folder <- "tests/testthat/data_metrictest"
files <- list.files(data_folder)
trajectory_data <- lapply(files, read_trajectory_ods)
trajectoryplot_data <- lapply(trajectory_data, plotLearnerData.ti.default)

all_trajs <- list(
  space_states = bind_rows(trajectoryplot_data %>% map(~ .$space_states)),
  space_lines = bind_rows(trajectoryplot_data %>% map(~ .$space_lines)),
  space_samples = bind_rows(trajectoryplot_data %>% map(~ .$space_samples))
)
class(all_trajs) <- "dyneval::ti_dimred_wrapper"

plotLearner.ti.default(all_trajs) + facet_wrap(~name)
plotLearner.ti.default(trajectoryplot_data[[4]])


traj <- trajectory_data[[4]]

task_emdist <- emdist(traj)
plot_emdist(traj, task_emdist)

# gr <- igraph::graph_from_data_frame(traj$state_network, directed = T, vertices = traj$state_names)
# milestone_distances <- igraph::distances(gr, weights = igraph::E(gr)$length)
# pct <- as.matrix(traj$state_percentages[,-1])
# rownames(pct) <- traj$state_percentages$id
#
#
#
#
# # can also implement myself https://people.cs.umass.edu/~mcgregor/papers/13-approx1.pdf
# cell_dists <- expand.grid(from = rownames(pct), to = rownames(pct)) %>%
#   mutate(dist = mapply(from, to, FUN = function(i, j) {
#     tr <- transport::transport(pct[i,], pct[j,], costm = milestone_distances, method = "revsimplex") %>%
#       mutate(dist = milestone_distances[cbind(from,to)], mult = mass * dist)
#     sum(tr$mult)
#   }))
# pheatmap::pheatmap(reshape2::acast(cell_dists, from~to, value.var = "dist"), cluster_rows = F, cluster_cols = F, annotation_col = as.data.frame(pct), annotation_row = as.data.frame(pct))
#
#
#
# x <- pct[1,]
# y <- pct[2,]
# tr <-  transport::transport(x, y, costm = milestone_distances, method = "shortsimplex")
# tr <-  transport::transport(x, y, costm = milestone_distances, method = "revsimplex")
# out <- tr %>%
#   mutate(dist = milestone_distances[cbind(from,to)], mult = mass * dist)
# out


# a <- c(100, 200, 80, 150, 50, 140, 170, 30, 10, 70)
# b <- c(60, 120, 150, 110, 40, 90, 160, 120, 70, 80)
# costm <- matrix(sample(1:20, 100, replace=TRUE), 10, 10)
# res <- transport::transport(a,b,costm) %>% mutate(dist = costm[cbind(from,to)], mult = mass * dist)
# distance <- sum(res$mult)

