% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/viz.R
\name{MarkerDimPlot}
\alias{MarkerDimPlot}
\title{Create a dimplot with markers}
\usage{
MarkerDimPlot(object, group_by, n = 20)
}
\arguments{
\item{object}{The Seurat object}

\item{group_by}{The column in the Seurat object meta.data to group data by (e.g. "Cell" or "seurat_clusters")}

\item{n}{The number of markers per identity to display}

\item{split_by}{The column in the Seurat object meta.data used to split each bar}
}
\value{
A plot object featuring a DimPlot next to a table of markers
}
\description{
Create a dimplot with markers
}
\examples{
MarkerDimPlot(seurat_obj, 'cell_type')
}
