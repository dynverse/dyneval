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

read_ods_trajectory <- function(filename) {
  structure_coltypes <- cols_only(from = col_character(), to = col_character(), length = col_double())

  structure <- my_read_ods(
    filename,
    sheet = "structure",
    col_types = structure_coltypes
  ) %>%
    mutate(
      from = paste0("milestone_", from),
      to = paste0("milestone_", to)
    )

  milestones <- sort(unique(c(structure$from, structure$to)))

  cells <- my_read_ods(
    filename,
    sheet = "cells",
    col_types = cols(.default = col_double())
  ) %>%
    select(one_of(milestones))

  for (milestone in milestones) {
    cells[,milestone] <- ifelse(is.na(cells[,milestone]), 0, cells[,milestone])
  }
  rownames(cells) <- paste0("cell_", seq_len(nrow(cells)))

  list(structure = structure, cells = cells)
}
