context("Getting all descriptions")

test_that("Descriptions can be retrieved", {
  tib <- get_descriptions()
  expect_that(tib, is_a("tbl"))

  lis <- get_descriptions(as_tibble = F)
  expect_that(lis, is_a("list"))

  for (descr in lis) {
    test_that(paste0("Description ", descr$name), {
      expect_lte(nchar(descr$short_name), 8)
    })
  }
})

test_that("Checking for dependencies does not produce an error", {
  expect_error(check_dependencies(), NA)
})


# run_method()
