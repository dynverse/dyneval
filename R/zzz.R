.onLoad <- function(libname, pkgname){
  cat("Loading discreteNameToValue in your global environment -- this is a dirty fix.\n")
  assign(
    "discreteNameToValue",
    ParamHelpers::discreteNameToValue,
    envir = .GlobalEnv
  )
}
