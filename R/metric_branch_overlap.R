
#' Compares the overlap between branches
#'
#' @param dataset A dataset
#' @param prediction A predicted model
#'
#' @importFrom dynfeature calculate_overall_feature_importance
calculate_branch_overlap <- function(dataset, prediction) {
  if (is.null(prediction) || is.null(dataset)) {
    lst(
      recovery_branches = 0,
      relevance_branches = 0,
      F1_branches = 0
    )
  } else {
    dataset <- dynwrap::simplify_trajectory(dataset, allow_self_loops = TRUE)
    prediction <- dynwrap::simplify_trajectory(prediction, allow_self_loops = TRUE)

    groups_dataset <- dataset %>% dynwrap::group_onto_trajectory_edges() %>% factor() %>% enframe("cell_id", "group_dataset")
    groups_prediction <- prediction %>% dynwrap::group_onto_trajectory_edges() %>% factor() %>% enframe("cell_id", "group_prediction")

    groups <- full_join(groups_dataset, groups_prediction, "cell_id")

    # calculate the size of the intersections and of each group separately
    intersections <- groups %>%
      filter(!is.na(group_dataset) & !is.na(group_prediction)) %>%
      group_by(group_dataset, group_prediction) %>%
      summarise(intersection = n()) %>%
      ungroup() %>%
      complete(group_dataset, group_prediction, fill = list(intersection = 0))

    n_dataset <- groups %>%
      filter(!is.na(group_dataset)) %>%
      group_by(group_dataset) %>%
      summarise(n_dataset = n())

    n_prediction <- groups %>%
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

    relevances <- jaccards %>%
      group_by(group_prediction) %>%
      arrange(-jaccard) %>%
      slice(1) %>%
      ungroup()

    # calculate the final scores
    zero_if_na <- function(x) if(is.na(x)) {0} else {x}

    lst(
      recovery_branches = mean(recoveries$jaccard) %>% zero_if_na(),
      relevance_branches = mean(relevances$jaccard) %>% zero_if_na(),
      F1_branches = calculate_harmonic_mean(recovery_branches, relevance_branches)
    )
  }
}
