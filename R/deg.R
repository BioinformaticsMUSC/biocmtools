#' Find Markers for all celltypes + conditions using Seurat FindMarkers
#'
#' @param seurat_obj The Seurat object
#' @param label_col The label (celltype, cluster, etc) column for organizing groups
#' @param condition_col The condition column (treatment, etc) used for comparing between
#' @param background_group The background group for use in the tests (usually "control" or similar)
#' @param assay The assay to use, should be 'RNA' or 'SCT'
#' @param test.use The type of test to use (e.g. "wilcox", "MAST", etc)
#' @param latent.vars Any variables to regress out
#' @return Table with all differentially expressed genes
#' @importFrom Seurat Idents FindMarkers 
#' @importFrom tidyr unite
#' @examples
#' FindAllDEGs(seurat_obj, label_col = "celltype", condition_group = "condition", background_group = "control", test.use = "wilcox")
#' @export
FindAllDEGs <- function (seurat_obj, 
                             label_col,
                             condition_col,
                             background_group,
                             test.use = "wilcox",
                             assay = "RNA",
                             latent.vars = NULL,
                         ...) {
  
  #make sure latent.vars is only used for specific tests
  if (!(test.use) %in% c("LR", "negbinom", "poisson", "MAST")) {
    if (!is.null(latent.vars)) {
      cat(paste0("Latent.vars cannot be used for test type: '", test.use, "'. Running without latent.vars...\n"))
    }
  }
  #create copy of seurat object
  tmp = seurat_obj
  
  #find and make sure there are only 2 conditions
  conditions = unique(tmp@meta.data[[condition_col]])
  if (length(conditions) != 2) {
    stop("You must have exactly 2 conditions in the 'condition_col'")
  }
  if (!(background_group) %in% conditions){
    stop(paste0("Background group: '", background_group,"' was not found in the condition_col."))
  }
  test_group = conditions[which(conditions != background_group)]
  #create label_condition column
  joint_column_name = paste(label_col, condition_col, sep = "_")
  tmp@meta.data <- tmp@meta.data |>
    tidyr::unite(!!as.name(joint_column_name), !!as.name(label_col), !!as.name(condition_col), remove=FALSE)
  
  all_markers = NULL
  for (label in unique(tmp@meta.data[[label_col]])) {
    cat("working on", label, "\n")
    Seurat::Idents(tmp) <- joint_column_name
    sub_markers <- Seurat::FindMarkers(tmp,
                               ident.1 = paste(label, test_group, sep = "_"),
                               ident.2 = paste(label, background_group, sep = "_"),
                               only.pos = FALSE, 
                               min.pct = 0.10, 
                               test.use = test.use,
                               logfc.threshold = 0.10, 
                               assay=assay,
                               latent.vars = latent.vars)
    sub_markers$diff <- sub_markers$pct.1 - sub_markers$pct.2
    sub_markers$celltype <-label
    
    all_markers <- rbind(all_markers, sub_markers)
    
  }
  
  return (all_markers)
}
