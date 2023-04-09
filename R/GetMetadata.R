#' Cluster-level cell type prediction
#'
#' @param rna A data.frame or matrix with genes listed as row names
#' @param ref.data Reference gene expression. If not specified, reference gene
#' expression data will be downloaded from the DISCO database.
#' @param ref.deg DEG list.
#' @param atlas Specify the atlas. If not specified, CELLiD will use all data.
#' @param n.predict Number of predicted cell type. Maximum value is 3
#' @param ref.path Specify the file path to save the downloaded reference data.
#' @param ncores Number of core used for CELLiD
#' @export
#' @importFrom jsonlite fromJSON
GetSampleMetadata <- function(project = NULL, tissue = NULL, disease = NULL, sample.type = NULL, cell.type = NULL) {

  tryCatch({
    message("Retrieving metadata from DISCO database")
    meta = fromJSON("https://immunesinglecell.org/toolkitapi/getSampleCtInfo")
  }, error = function(e){
    stop("Failed to retrieve metadata. Please try again. If the issue persists, please contact us at li_mengwei@immunol.a-star.edu.sg for assistance.")
  })

  message("Filtering metadata")

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
    tryCatch({
      message("Retrieving ontology from DISCO database")
      cell.ontology = fromJSON("https://immunesinglecell.org/toolkitapi/getCellOntology")
    }, error = function(e){
      stop("Failed to retrieve ontology Please try again. If the issue persists, please contact us at li_mengwei@immunol.a-star.edu.sg for assistance.")
    })
  }
  cell.type = cell.ontology$cell_name
  cell.type = grep(term, cell.type, ignore.case = T)
  return(cell.type)
}
