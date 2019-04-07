#' Compares the mapping between milestones or branches
#'
#' @param dataset A dataset
#' @param prediction A predicted model
#' @param grouping How to group the cells, either branches or milestones
#' @param simplify Whether to simplify the trajectory (allowing self loops)
#'
#' @keywords metric
#'
#' @importFrom dynfeature calculate_overall_feature_importance
calculate_mapping <- function(dataset, prediction, grouping = c("branches", "milestones"), simplify = TRUE) {
  grouping <- match.arg(grouping)

  if (is.null(prediction) || is.null(dataset)) {
    lst(
      recovery = 0,
      relevance = 0,
      F1 = 0
    )
  } else {
    if (simplify) {
      dataset <- dynwrap::simplify_trajectory(dataset, allow_self_loops = TRUE)
      prediction <- dynwrap::simplify_trajectory(prediction, allow_self_loops = TRUE)
    }

    if (grouping == "branches") {
      groups_dataset <- dataset %>% dynwrap::group_onto_trajectory_edges()
      groups_prediction <- prediction %>% dynwrap::group_onto_trajectory_edges()
    } else if (grouping == "milestones") {
      groups_dataset <- dataset %>% dynwrap::group_onto_nearest_milestones()
      groups_prediction <- prediction %>% dynwrap::group_onto_nearest_milestones()
    }

    groups_dataset <- groups_dataset %>% as.character() %>% enframe("cell_id", "group_dataset")
    groups_dataset_levels <- unique(groups_dataset$group_dataset) %>% na.omit()
    groups_prediction <- groups_prediction %>% as.character() %>% enframe("cell_id", "group_prediction")
    groups_prediction_levels <- unique(groups_prediction$group_prediction) %>% na.omit()

    groups <- full_join(groups_dataset, groups_prediction, "cell_id")

    # calculate the size of the intersections and of each group separately
    intersections <-
      groups %>%
      filter(!is.na(group_dataset), !is.na(group_prediction)) %>%
      group_by(group_dataset, group_prediction) %>%
      summarise(intersection = n()) %>%
      ungroup() %>%
      mutate(
        group_dataset = factor(group_dataset, levels = groups_dataset_levels),
        group_prediction = factor(group_prediction, levels = groups_prediction_levels)
      ) %>%
      complete(
        group_dataset, group_prediction,
        fill = list(intersection = 0)
      ) %>%
      mutate_if(is.factor, as.character)

    n_dataset <-
      groups %>%
      filter(!is.na(group_dataset)) %>%
      group_by(group_dataset) %>%
      summarise(n_dataset = n())

    n_prediction <-
      groups %>%
      filter(!is.na(group_prediction)) %>%
      group_by(group_prediction) %>%
      summarise(n_prediction = n())

    # now join and calculate the jaccard
    jaccards <- intersections %>%
      left_join(n_dataset, "group_dataset") %>%
      left_join(n_prediction, "group_prediction") %>%
      mutate(jaccard = intersection / (n_dataset + n_prediction - intersection))

    # calculate the recovery and relevance
    recoveries <- jaccards %>%
      group_by(group_dataset) %>%
      arrange(-jaccard) %>%
      slice(1) %>%
      ungroup()

    relevances <-
      jaccards %>%
      group_by(group_prediction) %>%
      arrange(-jaccard) %>%
      slice(1) %>%
      ungroup()

    # calculate the final scores
    zero_if_na <- function(x) if(is.na(x)) {0} else {x}

    lst(
      recovery = mean(recoveries$jaccard) %>% zero_if_na(),
      relevance = mean(relevances$jaccard) %>% zero_if_na(),
      F1 = calculate_harmonic_mean(recovery, relevance)
    )
  }
}

#' @rdname calculate_mapping
calculate_mapping_milestones <- function(dataset, prediction, simplify = TRUE) {
  mapping <- calculate_mapping(dataset, prediction, "milestones", simplify = simplify)
  names(mapping) <- paste0(names(mapping), "_milestones")
  mapping
}


#' @rdname calculate_mapping
calculate_mapping_branches <- function(dataset, prediction, simplify = TRUE) {
  mapping <- calculate_mapping(dataset, prediction, grouping = "branches", simplify = simplify)
  names(mapping) <- paste0(names(mapping), "_branches")
  mapping
}
