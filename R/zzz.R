.onAttach <- function(libname, pkgname) {

  tryCatch({
    if(packageVersion("DISCOtoolkit") == fromJSON("https://www.immunesinglecell.org/api/vishuo/getToolkitVersion")$version) {
      packageStartupMessage("Welcome to DISCOtoolkit.")
    } else {
      packageStartupMessage("There is a new version of DISCOtoolkit available. Please update to ensure that the functions work normally.")
    }
  }, error = function(e){
    packageStartupMessage("Welcome to DISCOtoolkit.")
  })
}

.onLoad <- function(libname, pkgname) {
  op <- options()

  op.disco <- list()
  tryCatch({
    op.disco <- list(
      disco.url = fromJSON("https://www.immunesinglecell.org/api/vishuo/getToolkitUrl")$url
    )
  }, error = function(e){
    op.disco <- list(
      disco.url = "https://disco.bii.a-star.edu.sg/"
    )
  })

  toset <- !(names(op.disco) %in% names(op))
  if (any(toset)) options(op.disco[toset])

  invisible()
}
