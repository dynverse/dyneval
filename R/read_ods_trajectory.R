#' read data from ods files
#'
#' Bugfix wrapper of the \code{\link[readODS]{read_ods}} function
#'
#' @import readODS
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

#' Read an ODS file containing a trajectory. TODO: explain format. Ask Robrecht about the format for now.
#'
#' @param filename the location of the ODS
#'
#' @import dplyr
#' @import tidyr
#' @import readr
#'
#' @export
read_ods_trajectory <- function(filename) {
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
    mutate(cellid = paste0("cell_", cellid)) %>%
    select(cellid, one_of(milestones))

  # fill in NAs
  for (milestone in milestones) {
    cells[,milestone] <- ifelse(is.na(cells[,milestone]), 0, cells[,milestone])
  }

  list(structure = structure, cells = cells)
}
