context("Getting all descriptions")

tib <- get_descriptions()
expect_that(tib, is_a("tbl"))

lis <- get_descriptions(as_tibble = F)
expect_that(lis, is_a("list"))

for (descr in lis) {
  testthat(paste0("Description ", descr$name), {
    expect_lte(str_length(descr$short_name), 8)
  })
}

check_dependencies()

# run_method()
