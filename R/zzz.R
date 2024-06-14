.onAttach <- function(libname, pkgname) {

  tryCatch({
    if(packageVersion("DISCOtoolkit") == fromJSON("http://www.immunesinglecell.org/api/vishuo/getToolkitVersion")$version) {
      packageStartupMessage("Welcome to DISCOtoolkit.")
    } else {
      packageStartupMessage("A new version of DISCOtoolkit is available. Please update to ensure all functions work properly.")
    }
  }, error = function(e){
    packageStartupMessage("Welcome to DISCOtoolkit.")
  })
}

.onLoad <- function(libname, pkgname) {
  op <- options()

  op.disco <-  list(
    disco_url = "https://immunesinglecell.org/disco_v3_api/",
    timeout = 6000
  )

  # tryCatch({
  #   op.disco <- list(
  #     disco_url = fromJSON("http://www.immunesinglecell.org/api/vishuo/getToolkitUrl")$url,
  #     timeout = 6000
  #   )
  # }, error = function(e){
  #   message("Fail to get url prefix, use default value")
  # })

  # toset <- !(names(op.disco) %in% names(op))
  # if (any(toset)) options(op.disco[toset])
  options(op.disco)
  invisible()
}


GetJson <- function(url, info.msg, error.msg){
  tryCatch({
    message(info.msg)
    return(fromJSON(url))
  }, error = function(e){
    stop(error.msg)
  })
}

