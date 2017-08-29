context("Plot functions")

# datasets <- source_file("helper-datasets-load.R", chdir = F)
datasets <- readRDS(paste0(tempdir(), "/dyneval_test_datasets.rds"))

test_that("perform dimred on trajectory", {
  for (dataseti in seq_len(nrow(datasets))) {
    dataset <- extract_row_to_list(datasets, dataseti)
    dimred <- plotdata_default(dataset, insert_phantom_edges = T)

    # check dimred$space_milestones
    milestone_ids <- dataset$milestone_ids
    space_milestones <- dimred$space_milestones
    expect_equal( nrow(space_milestones), length(milestone_ids) )
    expect_equal( space_milestones$id, milestone_ids )
    expect_true( all(c("id", "Comp1", "Comp2", "colour") %in% colnames(space_milestones)) )
    expect_true( all(is.finite(space_milestones$Comp1)) )
    expect_true( all(is.finite(space_milestones$Comp2)) )

    # check dimred$milestone_network
    milestone_network <- dataset$milestone_network
    space_lines <- dimred$space_lines
    expect_equal( nrow(space_lines), nrow(milestone_network) )
    expect_true( all(is.finite(space_lines$from.Comp1)) )
    expect_true( all(is.finite(space_lines$from.Comp2)) )
    expect_true( all(is.finite(space_lines$to.Comp1)) )
    expect_true( all(is.finite(space_lines$to.Comp2)) )
    expect_equal( sort(unique(c(space_lines$from, space_lines$to))), sort(milestone_ids) )

    # check dimred$space_samples
    space_samples <- dimred$space_samples
    cell_ids <- dataset$cell_ids
    expect_equal( nrow(space_samples), length(cell_ids) )

    # try different param
    dimred2 <- plotdata_default(dataset, insert_phantom_edges = F)
    # dimred2
  }
})
