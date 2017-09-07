context("Getting all descriptions")

tib <- get_descriptions()
expect_that(tib, is_a("tbl"))

lis <- get_descriptions(as_tibble = F)
expect_that(lis, is_a("list"))

for (descr in lis) {
  test_that(paste0("Description ", descr$name), {
    expect_lte(nchar(descr$short_name), 8)
  })
}

check_dependencies()

# run_method()
