.onLoad <- function(libname, pkgname) {
  # Assign estimate_theta to a variable in your package's namespace
  if (requireNamespace("sceptre", quietly = TRUE)) {
    assign("estimate_theta", sceptre:::estimate_theta, envir = parent.env(environment()))
  }
}
