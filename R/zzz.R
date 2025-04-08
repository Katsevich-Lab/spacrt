.onLoad <- function(libname, pkgname) {
  # Dynamically retrieve the unexported function from the sceptre namespace
  if (requireNamespace("sceptre", quietly = TRUE)) {
    estimate_theta <- get("estimate_theta", envir = asNamespace("sceptre"))
    assign("estimate_theta", estimate_theta, envir = parent.env(environment()))
  }
}
