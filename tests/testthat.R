library(testthat)
library(dyneval)
library(dynmethods)
library(dyntoy)
library(dynutils)
library(dplyr)
library(purrr)
library(mlrMBO)

Sys.setenv("R_TESTS" = "")

test_check("dyneval")

