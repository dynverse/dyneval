context("Score metrics")

test_that(paste0("Check lies network score"), {
  net1 <- tibble(from=c(1), to=c(1), directed=TRUE, length=1)
  net2 <- tibble(from=c(1), to=c(1), directed=TRUE, length=1)

  expect_equal(dyneval:::calculate_lies_network_score(net1, net2), 1)

  net1 <- tibble(from=c(1, 2, 2), to=c(2, 1, 1), directed=TRUE, length=1)
  net2 <- tibble(from=c(1, 2, 2), to=c(2, 3, 4), directed=TRUE, length=1)

  expect_equal(dyneval:::calculate_lies_network_score(net1, net2), 0.5)

  net1 <- tibble(from=c(1), to=c(1), directed=TRUE, length=1)
  net2 <- tibble(from=c(1), to=c(2), directed=TRUE, length=1)

  expect_lt(dyneval:::calculate_lies_network_score(net1, net2), 1)
})
