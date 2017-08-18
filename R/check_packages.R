#' @export
check_dependencies <- function() {
  functions <- lsf.str("package:dyneval")
  description_functions <- functions[grep("description_", functions)]
  for (descr_fun in description_functions) {
    descr <- do.call(descr_fun, list())
    required_packages <- c(descr$package_load, descr$package_installed)
    installed <- required_packages %in% rownames(installed.packages())
    if (any(!installed)) {
      warning(sQuote(descr$name), " requires the following packages still to be installed: ", paste(sQuote(required_packages[!installed]), collapse = ", "))
    }
  }
}
