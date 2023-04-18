.onAttach <- function(libname, pkgname) {

  tryCatch({
    if(packageVersion("DISCOtoolkit") == fromJSON("http://www.immunesinglecell.org/api/vishuo/getToolkitVersion")$version) {
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

  op.disco <-  list(
    disco.url = "https://www.immunesinglecell.org/toolkitapi/",
    timeout = 6000
  )

  tryCatch({
    op.disco <- list(
      disco.url = fromJSON("https://www.immunesinglecell.org/api/vishuo/getToolkitUrl")$url,
      timeout = 6000
    )
  }, error = function(e){
    message("Fail to get url prefix, use default value")
  })

  # toset <- !(names(op.disco) %in% names(op))
  # if (any(toset)) options(op.disco[toset])
  options(op.disco)
  invisible()
}
