context("Score metrics")

library(purrr)

test_that(paste0("Edge flip returns relevant results"), {
  topologies <- eval(formals(dyntoy:::generate_milestone_network)$model)

  milestone_networks <- list(
    linear1 = tribble(
      ~from, ~to, ~length, ~directed,
      "a",   "b", 1,       TRUE
    ),

    linear2 = tribble(
      ~from, ~to, ~length, ~directed,
      "a",   "b", 0.1,      FALSE,
      "b",   "c", 0.9,      FALSE
    ),

    linear3 = data_frame(
      from = letters %>% head(-1),
      to = letters %>% tail(-1),
      length = runif(length(from), 10, 30),
      directed = FALSE
    ),

    bifurcating1 = tribble(
      ~from, ~to, ~length, ~directed,
      "a",   "b", 1, TRUE,
      "b",   "c", 2, TRUE,
      "b",   "d", 3, TRUE
    ),

    bifurcating2 = tribble(
      ~from, ~to, ~length, ~directed,
      "b",   "a", 1, FALSE,
      "b",   "c", 2, FALSE,
      "b",   "d", 3, FALSE,
      "a",   "x", 4, FALSE,
      "c",   "y", 5, FALSE,
      "d",   "z", 6, FALSE
    ),

    cycle1 = tribble(
      ~from, ~to, ~length, ~directed,
      "a", "b", 1, TRUE,
      "b", "c", 2, TRUE,
      "c", "a", 3, TRUE
    ),

    cycle2 = data_frame(
      from = letters,
      to = c(letters[-1], letters[[1]]),
      length = runif(length(letters), 4, 10),
      directed = FALSE
    )
  )

  networks <- tibble(
    name = names(milestone_networks),
    milestone_network = milestone_networks,
    topology = gsub("\\d*$", "", name)
  )

  design <- crossing(networks, networks)

  scores <- design %>% as.list() %>% pmap(function(milestone_network, milestone_network1, ...) {
    score <- calculate_edge_flip(milestone_network, milestone_network1)
    tibble(score = score)
  }) %>% bind_rows()
  design <- bind_cols(design, scores)

  expect_type(design$score, "double")

  same_topology_scores <- design %>% mutate(test = (topology == topology1) == (abs(score - 1) < 1e-8)) %>% pull(test)

  expect_true(all(same_topology_scores))
})




test_that(paste0("Correlation returns relevant results"), {
  # test for when cells are on one point -> correlation == 0
  cell_ids <- letters[1:10]
  dataset <- dynwrap::wrap_data(cell_ids = cell_ids) %>%
    dynwrap::add_trajectory(
      milestone_network = tribble(
        ~from, ~to, ~length, ~directed,
        "A", "B", 1, TRUE
      ),
      progressions = tibble(
        cell_id = cell_ids,
        from = "A",
        to = "B",
        percentage = 0
      )
    ) %>%
    dynwrap::add_cell_waypoints()

  scores <- calculate_metrics(dataset = dataset, model = dataset, metrics = "correlation")

  expect_equal(scores$correlation, 0)
})
