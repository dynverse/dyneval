library(testthat)
library(dyneval)
library(dynwrap)
library(dyntoy)
library(dynutils)
library(dplyr)
library(purrr)
library(tidyr)

Sys.setenv("R_TESTS" = "")

test_check("dyneval")

