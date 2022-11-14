#' Assess batch correction by Harmony by printing out side by side dim plots before and after correction
#'
#' @param object The Seurat object
#' @param replicate_col The column in the Seurat object that distinguishes replicates (eg. "patient")
#' @return plot
#' @importFrom Seurat DimPlot
#' @importFrom ggplot2 ggtitle
#' @import patchwork 
#' @examples
#' CompareBatchEffects(seurat_obj, 'patient')
#' @export
CompareBatchEffects <- function(object, batch_col) {
  if (("pca" %in% names(object@reductions)) & ("harmony" %in% names(object@reductions))) {
    p1 <- DimPlot(object, reduction = 'pca', group.by = batch_col) + ggtitle("Before Batch Effect Correction")
    p2 <- DimPlot(object, reduction = 'harmony', group.by = batch_col) + ggtitle("After Batch Effect Correction")
    return (p1 + p2)
  }
  
}