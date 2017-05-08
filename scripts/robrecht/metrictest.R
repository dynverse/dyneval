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
  space_states = trajectoryplot_data %>% map_df(~ .$space_states),
  space_lines = trajectoryplot_data %>% map_df(~ .$space_lines),
  space_samples = trajectoryplot_data %>% map_df(~ .$space_samples)
)
class(all_trajs) <- "dyneval::ti_dimred_wrapper"

plotLearner.ti.default(all_trajs) + facet_wrap(~name)

out_folder <- "scratch/output_metrictest"
dir.create(out_folder)

pngsz <- 1200
for (ix in seq_along(trajectory_data)) {
  cat("Processing test ", ix, "\n", sep="")

  file_traj <- paste0(out_folder, "/test_", ix, "_traj.png")
  file_hm <- paste0(out_folder, "/test_", ix, "_heatmap.png")
  file_comb <- paste0(out_folder, "/test_", ix, "_comb.png")

  traj <- trajectory_data[[ix]]
  plotdata <- trajectoryplot_data[[ix]]

  png(file_traj, pngsz, pngsz, res = 300)
  print(plotLearner.ti.default(plotdata))
  dev.off()

  emdist <- compute_emlike_dist(traj)

  plot_emdist(traj, emdist,
              filename = file_hm,
              width = pngsz/300, height = pngsz/300)

  system(paste0("convert ", file_traj, " ", file_hm, " -append ", file_comb))
  file.remove(file_traj)
  file.remove(file_hm)
}

in_files <- paste(paste0(out_folder, "/test_", seq_along(trajectory_data), "_comb.png"), collapse = " ")

system(paste0("convert ", in_files, " +append ", out_folder, "/test_comb.png"))

