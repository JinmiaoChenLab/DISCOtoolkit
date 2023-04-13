
#' Get metadata from DISCO and apply filters on it
#'
#' @param sample sample ID(s) of interest, can be a string or vector
#' @param project project ID(s) of interest, can be a string or vector
#' @param tissue tissue(s) of interest, can be a string or vector
#' @param disease disease(s) of interest, can be a string or vector
#' @param platform platform(s) of interest, can be a string or vector
#' @param sample.type sample type(s) of interest, can be a string or vector
#' @param cell.type cell type(s) of interest, can be a string or vector
#' @param cell.type.confidence  Filter the results based on the confidence of cell type prediction, which can be categorized as high, medium, or all
#' @param include.cell.type.children Whether the cell type is an exact match or includes its child elements
#' @param min.cell.per.sample Only samples with a cell count greater than the minimum cell count per sample will be retained.
#' @export
#' @importFrom jsonlite fromJSON
FilterDiscoMetadata <- function(sample = NULL,
                                project = NULL,
                                tissue = NULL,
                                disease = NULL,
                                platform = NULL,
                                sample.type = NULL,
                                cell.type = NULL,
                                cell.type.confidence = "medium",
                                include.cell.type.children = T,
                                min.cell.per.sample = 100) {

  filter.data = list(
    sample.metadata = NULL,
    cell.type.metadata = NULL,
    sample.count = NULL,
    cell.count = NULL,
    filter =  list(
      sample = sample, project = project, tissue = tissue, disease = disease, platform = platform,
      sample.type = sample.type, cell.type = cell.type, cell.type.confidence = cell.type.confidence,
      include.cell.type.children = include.cell.type.children, min.cell.per.sample = min.cell.per.sample
    )
  )

  metadata = GetDiscoMetadata()

  message("Filtering sample")

  if (is.null(sample) == F) {
    metadata = metadata[which(metadata$sampleId %in% sample),]
  }

  if (is.null(project) == F) {
    metadata = metadata[which(metadata$projectId %in% project),]
  }

  if (is.null(tissue) == F) {
    metadata = metadata[which(metadata$tissue %in% tissue),]
  }

  if (is.null(platform) == F) {
    metadata = metadata[which(metadata$platform %in% platform),]
  }

  if (is.null(disease) == F) {
    metadata = metadata[which(metadata$disease %in% disease),]
  }

  if (is.null(sample.type) == F) {
    metadata = metadata[which(metadata$sampleType %in% sample.type),]
  }

  if (nrow(metadata) == 0) {
    stop("Sorry, no samples passed the applied filters.")
  }

  sample.ct.info = GetSampleCtInfo()


  retain.field = c("sampleId", "projectId", "sampleType", "anatomicalSite", "disease",
                   "tissue", "platform", "ageGroup", "age", "gender", "cellSorting",
                   "diseaseSubtype", "diseaseStage", "treatment", "md5")
  metadata = metadata[,retain.field]

  if (is.null(cell.type)) {

    sample.ct.info = sample.ct.info[which(sample.ct.info$sampleId %in% metadata$sampleId),]
    sample.cell.count = aggregate(sample.ct.info$cellNumber, by=list(sample=sample.ct.info$sampleId), FUN=sum)
    rownames(sample.cell.count) = sample.cell.count$sample
    metadata$cell.number = sample.cell.count[rownames(metadata),"x"]

    metadata = metadata[which(metadata$cell.number > min.cell.per.sample),]
    if (nrow(metadata) == 0) {
      stop("Sorry, no samples passed the applied filters.")
    }
    sample.ct.info = sample.ct.info[which(sample.ct.info$sampleId %in% metadata$sampleId),]

    filter.data[["sample.metadata"]] = metadata
    filter.data[["cell.type.metadata"]] = sample.ct.info
    filter.data[["sample.count"]] = nrow(metadata)
    filter.data[["cell.count"]] = sum(metadata$cell.number)
    message(paste0(filter.data[["sample.count"]], " samples and ", filter.data[["cell.count"]]," cells were found"))
    return(filter.data)
  }

  if (include.cell.type.children) {
    cell.type = GetCellTypeChildren(cell.type)
  }

  sample.ct.info = sample.ct.info[which(sample.ct.info$cellType %in% cell.type),]
  sample.ct.info = sample.ct.info[which(sample.ct.info$sampleId %in% metadata$sampleId),]

  if (!(cell.type.confidence %in% c("high", "medium", "all"))) {
    stop("cell.type.confidence can only be high, medium, or all")
  }

  if (cell.type.confidence == "high") {
    sample.ct.info = sample.ct.info[which(sample.ct.info$cellTypeScore >= 0.8),]
  } else if  (cell.type.confidence == "medium") {
    sample.ct.info = sample.ct.info[which(sample.ct.info$cellTypeScore >= 0.6),]
  }


  if (nrow(sample.ct.info) == 0) {
    stop("Sorry, no samples passed the applied filters.")
  }

  metadata = metadata[which(metadata$sampleId %in% sample.ct.info$sampleId),]
  sample.cell.count = aggregate(sample.ct.info$cellNumber, by=list(sample=sample.ct.info$sampleId), FUN=sum)

  rownames(sample.cell.count) = sample.cell.count$sample
  metadata$cell.number = sample.cell.count[rownames(metadata),"x"]

  metadata = metadata[which(metadata$cell.number > min.cell.per.sample),]
  if (nrow(metadata) == 0) {
    stop("Sorry, no samples passed the applied filters.")
  }
  sample.ct.info = sample.ct.info[which(sample.ct.info$sampleId %in% metadata$sampleId),]

  filter.data[["sample.metadata"]] = metadata
  filter.data[["cell.type.metadata"]] = sample.ct.info
  filter.data[["sample.count"]] = nrow(metadata)
  filter.data[["cell.count"]] = sum(metadata$cell.number)

  message(paste0(filter.data[["sample.count"]], " samples and ", filter.data[["cell.count"]]," cells were found"))

  return(filter.data)
}

#' Check children of input cell type
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
#' @export
GetCellTypeChildren <- function(cell.type, cell.ontology = NULL) {
  if(is.null(cell.ontology)) {
    tryCatch({
      message("Retrieving ontology from DISCO database")
      cell.ontology = fromJSON("https://immunesinglecell.org/toolkitapi/getCellOntology")
    }, error = function(e){
      stop("Failed to retrieve ontology Please try again. If the issue persists, please contact us at li_mengwei@immunol.a-star.edu.sg for assistance.")
    })
  }

  children = c(cell.type)

  while(length(cell.type) > 0) {
    children = append(children, cell.ontology$cell_name[which(cell.ontology$parent %in% cell.type)])
    cell.type = cell.ontology$cell_name[which(cell.ontology$parent %in% cell.type)]
  }

  children = unique(children)
  return(children)
}

#' Get metadata of all samples in DISCO
#'
#' @export
GetDiscoMetadata <- function() {

  metadata = GetJson(
    url = "https://www.immunesinglecell.org/api/vishuo/sample/all",
    info.msg = "Retrieving metadata from DISCO database",
    error.msg = "Failed to retrieve metadata. Please try again. If the issue persists, please contact us at li_mengwei@immunol.a-star.edu.sg for assistance."
  )
  metadata = metadata[which(metadata$processStatus == "QC pass"),]
  rownames(metadata) = metadata$sampleId
  return(metadata)
}


#' Check children of input cell type
#'
#' @param term A partial or complete name of a cell type
#' @param cell.ontology The cell type ontology. If not specified, the cell type ontology will be retrieved from the DISCO database.
#' @return List of cell type
#' @examples
#' # Retrieve cell types whose names contain the term 'blast'
#' GetCellTypeChildren(cell.type = c("blast"))
#' @export
FindCellType <- function(term = "", cell.ontology = NULL) {
  if(is.null(cell.ontology)) {
    cell.ontology = GetJson(
      url = "https://immunesinglecell.org/toolkitapi/getCellOntology",
      info.msg = "Retrieving ontology from DISCO database",
      error.msg = "Failed to retrieve ontology Please try again. If the issue persists, please contact us at li_mengwei@immunol.a-star.edu.sg for assistance."
    )
  }
  cell.type = cell.ontology$cell_name
  cell.type = grep(term, cell.type, ignore.case = T, value = T)
  return(cell.type)
}


GetJson <- function(url, info.msg, error.msg){
  tryCatch({
    message(info.msg)
    return(fromJSON(url))
  }, error = function(e){
    stop(error.msg)
  })
}

#' Get cell type information of sample
GetSampleCtInfo <- function() {
  sample.ct.info = GetJson(
    url = "https://immunesinglecell.org/toolkitapi/getSampleCtInfo",
    info.msg = "Retrieving cell type information of each sample from DISCO database",
    error.msg = "Failed to retrieve cell type information. Please try again. If the issue persists, please contact us at li_mengwei@immunol.a-star.edu.sg for assistance."
  )
  return(sample.ct.info)
}


#' List all potential values of a metadata field
#'
#' @export
ListMetadataItem <- function(field){
  metadata = GetDiscoMetadata()
  if (field %in% colnames(metadata)) {
    return(unique(metadata[,field]))
  } else {
    stop(paste0("DISCO data don't contain ", field, " field"))
  }
}


