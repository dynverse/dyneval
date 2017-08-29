context("Data IO")

# datasets <- source_file("helper-datasets-load.R", chdir = F)
datasets <- readRDS(paste0(tempdir(), "/dyneval_test_datasets.rds"))

test_that("Loading datasets", {
  expect_that( is_tibble(datasets), is_true() )

  required_cols <- c("id", "cell_ids", "milestone_ids", "milestone_network", "milestone_percentages", "progressions", "counts", "geodesic_dist", "special_cells")
  expect_that( all(required_cols %in% colnames(datasets)), is_true() )
})


for (dataseti in seq_len(nrow(datasets))) {
  dataset <- extract_row_to_list(datasets, dataseti)

  test_that(paste0("Evaluating with ", dataset$id), {


    expect_true( is.character(dataset$id) )

    cell_ids <- dataset$cell_ids
    milestone_ids <- dataset$milestone_ids

    expect_true( is.character(cell_ids) )
    expect_true( is.character(milestone_ids) )

    milestone_network <- dataset$milestone_network
    expect_true( all(milestone_network$from %in% milestone_ids) )
    expect_true( all(milestone_network$to %in% milestone_ids) )
    expect_true( all(milestone_network$length > 0) )

    milestone_percentages <- dataset$milestone_percentages
    expect_equal( sort(unique(milestone_percentages$cell_id)), sort(cell_ids) )
    expect_equal( sort(unique(milestone_percentages$milestone_id)), sort(milestone_ids) )
    # todo: check whether a cell is in a tent

    mp_summ <- milestone_percentages %>% group_by(cell_id) %>% summarise(sum = sum(percentage))
    expect_true( all(mp_summ$sum == 1) )

    progressions <- dataset$progressions
    expect_equal( sort(unique(progressions$cell_id)), sort(cell_ids) )
    expect_equal( sort(unique(c(progressions$from, progressions$to))), sort(milestone_ids) )

    pr_summ <- progressions %>% group_by(cell_id) %>% summarise(sum = sum(percentage))
    expect_true( all(mp_summ$sum <= 1) )

    counts <- dataset$counts
    expect_true( is.numeric(counts) )
    expect_true( is.matrix(counts) )
    expect_equal( rownames(counts), cell_ids )

    geodesic_dist <- dataset$geodesic_dist
    expect_true( is.matrix(geodesic_dist) )
    expect_true( is.numeric(geodesic_dist) )
    expect_equal( rownames(geodesic_dist), cell_ids )
    expect_equal( colnames(geodesic_dist), cell_ids )
    expect_true( all(is.finite(geodesic_dist)) )
    expect_true( all(geodesic_dist >= 0) )

    special_cells <- dataset$special_cells
    expect_true( all(sapply(special_cells, function(x) all(x %in% cell_ids))) )
  })
}
