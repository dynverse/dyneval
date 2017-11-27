library(testthat)
library(dyneval)
library(dynmethods)
library(dyntoy)
library(dynutils)
library(dplyr)
library(purrr)
library(mlrMBO)
library(tidyr)

Sys.setenv("R_TESTS" = "")

test_check("dyneval")

