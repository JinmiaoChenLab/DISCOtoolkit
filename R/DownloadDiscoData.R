#' Download data from DISCO database
#'
#' @param cell.type The cell type of the input, which may be a string or a list of string
#' @param cell.ontology The cell type ontology. If not specified, the cell type ontology will be retrieved from the DISCO database.
#' @return List of cell type
#' @examples
#' # Get children of B cell
#' GetCellTypeChildren(cell.type = c("B cell"))
#'
#' # Get children of B cell and Macrophage
#' GetCellTypeChildren(cell.type = c("B cell", "Macrophage"))
#' @importFrom tools md5sum
#' @export
DownloadDiscoData <- function(metadata, output.dir = "DISCOtmp") {

  error.samples = c()

  tryCatch({
    if (!dir.exists(output.dir)) {
      dir.create(output.dir)
    }
  }, error = function(e){
    stop("The output directory cannot be created.")
  })


  if (is.null(metadata$filter$cell.type)) {
    samples = metadata$sample.metadata
    for (i in 1:nrow(samples)) {
      output.file = paste0(output.dir, "/", samples$sampleId[i], ".rds")
      if (file.exists(output.file) & md5sum(output.file) == samples$md5[i]) {
        message(paste0(samples$sampleId[i], " has been downloaded before. Ignore..."))
      } else {
        message(paste0("Downloading data of ", samples$sampleId[i]))
        tryCatch({
          download.file(
            url = paste0("https://immunesinglecell.org/toolkitapi/getRdsBySample?sample=", samples$sampleId[i]),
            destfile = output.file
          )
        }, error = function(e) {
          error.samples = append(error.samples, samples$sampleId[i])
        })
      }
    }
  } else {
    samples = metadata$cell.type.metadata
    samples$sampleId = paste0(samples$sampleId, "_", samples$cluster)
    for (i in 1:nrow(samples)) {
      output.file = paste0(output.dir, "/", samples$sampleId[i], ".rds")
      if (file.exists(output.file) & md5sum(output.file) == samples$md5[i]) {
        message(paste0(samples$sampleId[i], " has been downloaded before. Ignore..."))
      } else {
        message(paste0("Downloading data of ", samples$sampleId[i]))
        tryCatch({
          download.file(
            url = paste0("https://immunesinglecell.org/toolkitapi/getRdsBySampleCt?sample=", samples$sampleId[i]),
            destfile = output.file
          )
        }, error = function(e) {
          error.samples = append(error.samples, samples$sampleId[i])
        })
      }
    }
  }

  return(error.samples)
}

