library(future)

test_fun <- function(method_sleep_time, wait_time) {
  start_time <- Sys.time()

  fut <- future({
    Sys.sleep(method_sleep_time)
    paste0("Seconds slept: ", method_sleep_time)
  }, evaluator = plan("multisession"))

  while (!resolved(fut) && difftime(Sys.time(), start_time, units = "secs") < wait_time) {
    Sys.sleep(1)
  }

  if (resolved(fut)) {
    value(fut)
  } else {
    future:::ClusterRegistry("stop")
    NA
  }
}

test_fun(4, 10)
test_fun(8, 10)
test_fun(10, 10)
test_fun(15, 10)

seconds <- c(seq(1, 20, length.out = 8), rep(5, 40))
unlist(parallel::mclapply(seconds, mc.cores = 8, test_fun, wait_time = 10), recursive = F)
