context("Getting all descriptions")

tib <- get_descriptions()
expect_that(tib, is_a("tbl"))

lis <- get_descriptions(as_tibble = F)
expect_that(lis, is_a("list"))

check_dependencies()

# run_method()
