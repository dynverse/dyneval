context("Plot functions")

tasks <- dyneval::generate_toy_datasets()

for (taski in seq_len(nrow(tasks))) {
  task <- extract_row_to_list(tasks, taski)

  test_that(paste0("Perform dimred on trajectory with task ", task$id), {
    g <- plot_default(task)
    expect_is(g$layers[[1]], "ggproto")
  })
}
