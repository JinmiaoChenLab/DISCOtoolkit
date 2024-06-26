% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GetMetadata.R
\name{FindCellType}
\alias{FindCellType}
\title{Check children of input cell type}
\usage{
FindCellType(term = "", cell_ontology = NULL)
}
\arguments{
\item{term}{A partial or complete name of a cell type}

\item{cell_ontology}{The cell type ontology. If not specified, the cell type ontology will be retrieved from the DISCO database.}
}
\value{
List of cell type
}
\description{
Check children of input cell type
}
\examples{
# Retrieve cell types whose names contain the term 'blast'
GetCellTypeChildren(cell_type = c("blast"))
}
