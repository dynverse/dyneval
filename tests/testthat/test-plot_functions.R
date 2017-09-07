context("Plot functions")

tasks <- readRDS(paste0(tempdir(), "/dyneval_test_datasets.rds"))

for (taski in seq_len(nrow(tasks))) {
  task <- extract_row_to_list(tasks, taski)

  test_that(paste0("Perform dimred on trajectory with task ", task$id), {
    g <- plot_default(task)
    expect_is(g, "ggplot")

    pdf("/dev/null")
    print(g)
    dev.off()

    prediction <- task
    cell_id_map <- setNames(sample(prediction$cell_ids), prediction$cell_ids)
    prediction$milestone_percentages$cell_id <- cell_id_map[prediction$milestone_percentages$cell_id]
    prediction$progressions$cell_id <- cell_id_map[prediction$progressions$cell_id]

    g <- plot_combined(task, prediction)
    expect_is(g, "ggplot")

    pdf("/dev/null")
    print(g)
    dev.off()

    pdf("/dev/null")
    ph <- plot_emdist(task, task$geodesic_dist)
    dev.off()
    expect_is(ph$gtable, "gtable")
  })
}
