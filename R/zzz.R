.onLoad <- function(libname, pkgname){
  packageStartupMessage("Loading discreteNameToValue in your global environment -- this is a dirty fix.")
  suppressWarnings({
    requireNamespace("ParamHelpers")
    assign(
      "discreteNameToValue",
      ParamHelpers::discreteNameToValue,
      envir = .GlobalEnv
    )
  })
}
