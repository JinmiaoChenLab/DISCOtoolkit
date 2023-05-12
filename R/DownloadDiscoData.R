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
        error.samples = tryCatch({
          download.file(
            url = paste0(getOption("disco.url"), "/getRdsBySample?sample=", samples$sampleId[i], "&project=", samples$projectId[i]),
            destfile = output.file, method = "curl"
          )
          if (!(md5sum(output.file) == samples$md5[i])) {
            error.samples = append(error.samples, samples$sampleId[i])
            unlink(output.file)
            message(paste(samples$sampleId[i]), " md5 check failed")
          }
        }, error = function(e) {
          error.samples = append(error.samples, samples$sampleId[i])
          return(error.samples)
        })
      }
    }
    if (length(error.samples) > 0) {
      metadata$sample.metadata = metadata$sample.metadata[error.samples,,drop=F]
      metadata$cell.type.metadata = metadata$cell.type.metadata[which(metadata$cell.type.metadata$sampleId %in% error.samples),,drop=F]
      metadata$cell.count = sum(metadata$cell.type.metadata$cellNumber)
      metadata$sample.count = length(error.samples)
      return(metadata)
    }
  } else {
    samples = metadata$cell.type.metadata
    samples$sample = samples$sampleId
    samples$sampleId = paste0(samples$sampleId, "_", samples$cluster)
    for (i in 1:nrow(samples)) {
      output.file = paste0(output.dir, "/", samples$sampleId[i], ".rds")
      if (file.exists(output.file) & md5sum(output.file) == samples$md5[i]) {
        message(paste0(samples$sampleId[i], " has been downloaded before. Ignore..."))
      } else {
        message(paste0("Downloading data of ", samples$sampleId[i]))
        error.samples = tryCatch({
          download.file(
            url = paste0(getOption("disco.url"),"/getRdsBySampleCt?sample=", samples$sampleId[i], "&project=", metadata$sample.metadata[samples$sample[i], "projectId"]),
            destfile = output.file, method = "curl"
          )
          if (!(md5sum(output.file) == samples$md5[i])) {
            error.samples = append(error.samples, samples$sampleId[i])
            unlink(output.file)
            message(paste(samples$sampleId[i]), " md5 check failed")
          }
        }, error = function(e) {
          error.samples = append(error.samples, samples$sampleId[i])
          return(error.samples)
        })
      }
    }

    if (length(error.samples) > 0) {
      metadata$cell.type.metadata = metadata$cell.type.metadata[which(samples$sampleId %in% error.samples),,drop=F]
      metadata$sample.metadata = metadata$sample.metadata[which(metadata$sample.metadata$sampleId %in% unique(metadata$cell.type.metadata$sampleId)),,drop=F]

      metadata$cell.count = sum(metadata$cell.type.metadata$cellNumber)
      metadata$sample.count = length(unique(metadata$cell.type.metadata$sampleId))

      sample.cell.count = aggregate(metadata$cell.type.metadata$cellNumber, by=list(sample=metadata$cell.type.metadata$sampleId), FUN=sum)
      rownames(sample.cell.count) = sample.cell.count$sample
      metadata$sample.metadata$cell.number = sample.cell.count[rownames(metadata$sample.metadata),"x"]

      return(metadata)
    }
  }

  return(NULL)
}

