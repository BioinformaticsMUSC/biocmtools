#' Prepare Seurat Object for AnnData Conversion
#'
#' This function exports a Seurat object to multiple CSV and MTX files that can
#' be used to reconstruct the data as an AnnData object in Python. It extracts
#' metadata, count matrices, dimensionality reductions, and gene information.
#'
#' @param seurat A Seurat object to be prepared for AnnData conversion
#' @param output_dir Character string specifying the output directory. If NULL,
#'   uses the current working directory. Directory will be created if it doesn't exist.
#'
#' @details The function exports the following files:
#' \itemize{
#'   \item{metadata.csv: Cell metadata including barcodes and UMAP coordinates (if available)}
#'   \item{counts.mtx: Gene expression count matrix in Matrix Market format}
#'   \item{sct_counts.mtx: SCT-transformed counts (if SCT assay is present)}
#'   \item{pca.csv: PCA cell embeddings (if available)}
#'   \item{harmony.csv: Harmony cell embeddings (if available)}
#'   \item{gene_names.csv: Gene names corresponding to the count matrix rows}
#' }
#'
#' @return Invisible NULL. Function is called for its side effects (file writing).
#'
#' @examples
#' \dontrun{
#' # Prepare a Seurat object for AnnData conversion
#' prep_seurat_for_anndata(my_seurat, output_dir = "./anndata_files")
#' }
#'
#' @export
#' @importFrom Matrix writeMM
#' @importFrom utils write.csv write.table
prep_seurat_for_anndata <- function(seurat, output_dir = NULL) {
  
  # Set output directory
  if (is.null(output_dir)) {
    output_dir <- getwd()
  }
  if (!(dir.exists(output_dir))) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Store original working directory and set new one
  original_wd <- getwd()
  setwd(output_dir)
  on.exit(setwd(original_wd))
  
  cat('Preparing seurat object for anndata conversion...\n')
  
  # Add barcode information to metadata
  seurat$barcode <- colnames(seurat)
  
  # Add UMAP coordinates if available
  if ("umap" %in% names(seurat@reductions)) {
    seurat$UMAP_1 <- seurat@reductions$umap@cell.embeddings[,1]
    seurat$UMAP_2 <- seurat@reductions$umap@cell.embeddings[,2]
  }
  
  # Write metadata
  cat('Writing metadata.csv...\n')
  write.csv(seurat@meta.data, file='metadata.csv', quote=FALSE, row.names=FALSE)
  
  # Write gene counts matrix
  cat('Creating gene counts matrix counts.mtx...\n')
  counts_matrix <- Seurat::LayerData(seurat, assay="RNA", layer="counts")
  Matrix::writeMM(counts_matrix, file='counts.mtx')
  
  # Write SCT counts if available
  if ("SCT" %in% names(seurat@assays)) {
    cat("SCT assay detected. Creating SCT counts matrix...\n")
    sct_counts_matrix <- Seurat::LayerData(seurat, assay="SCT", layer="counts")
    Matrix::writeMM(sct_counts_matrix, file='sct_counts.mtx')
  }
  
  # Write PCA embeddings if available
  if ("pca" %in% names(seurat@reductions)) {
    cat('Writing pca.csv...\n')
    write.csv(seurat@reductions$pca@cell.embeddings, file='pca.csv', 
              quote=FALSE, row.names=FALSE)
  }
  
  # Write Harmony embeddings if available
  if ("harmony" %in% names(seurat@reductions)) {
    cat('Harmony reduction detected. Writing harmony.csv...\n')
    write.csv(seurat@reductions$harmony@cell.embeddings, file='harmony.csv', 
              quote=FALSE, row.names=FALSE)
  }

  # Write gene names
  cat('Writing gene names table...\n')
  write.table(
    data.frame('gene'=rownames(counts_matrix)), file='gene_names.csv',
    quote=FALSE, row.names=FALSE, col.names=FALSE
  )
  
  cat('Done!\n')
  cat(paste('Files for Anndata have been saved to:', output_dir, '\n'))
  
  invisible(NULL)
}