.onLoad <- function(libname, pkgname){
  requireNamespace("ParamHelpers")
  packageStartupMessage("Loading discreteNameToValue in your global environment -- this is a dirty fix.")
  assign(
    "discreteNameToValue",
    ParamHelpers::discreteNameToValue,
    envir = .GlobalEnv
  )
}
