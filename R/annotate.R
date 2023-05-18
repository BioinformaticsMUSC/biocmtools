# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/blob/master/LICENSE)
# Written by Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi>, June 2021

#' Run sc Type on a seurat object
#'
#' @param seurat The Seurat object
#' @param tissue The tissue within sc-type to use
#' @param idents The clustering column to use as input for sc-type
#' @param assay The assay to use, should be 'RNA' or 'SCT'
#' @param new_col_name The name of the new column in the seurat object with sc-type predictions
#' @return Seurat object
#' @import Seurat 
#' @import dplyr
#' @import HGNChelper
#' @examples
#' RunSCtype(seurat_obj, tissue = "Brain", idents = "seurat_clusters")
#' @export
RunSCtype <- apply_scType <- function (seurat, tissue, idents = NULL, assay = NULL, new_col_name = NULL, verbose = TRUE) {
  if (is.null(tissue)) {
    stop('must select a tissue to proceed')
  }
  if (is.null(idents)) {
    Idents(seurat) <- "seurat_clusters"
  } else {
    if (idents %in% colnames(seurat@meta.data)) {
      Idents(seurat) <- idents
    } else {
      stop('"seurat_clusters" column not in metadata - please provide a clustering column')
    }
  }
  if (is.null(assay)) {
    assay <- "RNA"
  }
  if (is.null(new_col_name)) {
    new_col_name = paste0('sctype_', tissue, '_prediction')
  }
  source(
    "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"
  )
  source(
    "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R"
  )
  
  gs_list = gene_sets_prepare(
    "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx",
    tissue
  )
  
  es.max = sctype_score(
    scRNAseqData = seurat[['RNA']]@scale.data,
    scaled = TRUE,
    gs = gs_list$gs_positive,
    gs2 = gs_list$gs_negative
  )
  
  cL_resutls = do.call("rbind", lapply(unique(seurat@active.ident), function(cl) {
    es.max.cl = sort(rowSums(es.max[, rownames(seurat@meta.data[seurat@active.ident ==
                                                                  cl,])]), decreasing = !0)
    head(data.frame(
      cluster = cl,
      type = names(es.max.cl),
      scores = es.max.cl,
      ncells = sum(seurat@active.ident == cl)
    ),
    10)
  }))
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
  
  #if there are multipe cell assignments to 1 column
  sctype_scores <- distinct(sctype_scores, cluster, .keep_all = TRUE)
  
  #Map to new column
  current.cluster.ids <- sctype_scores$cluster
  new.cluster.ids <- as.character(sctype_scores$type)
  names(new.cluster.ids) <- current.cluster.ids
  seurat[[new_col_name]] <- recode(seurat@active.ident, !!!new.cluster.ids)
  
  if (isTRUE(verbose)) {
    print(sctype_scores)
  }
  
  return (seurat)
}

#' Assess batch correction by Harmony by printing out side by side dim plots before and after correction
#'
#' @param seurat The Seurat object
#' @param cluster_col Clusters in seurat object to use, default is 'seurat_clusters'
#' @param cell_type_col Cell type column in seurat, default is "sctype_id"
#' @return data.frame
#' @import Seurat 
#' @examples
#' CreateLabelFile(seurat)
#' @export
CreateLabelFile <- function (seurat, cluster_col = 'seurat_clusters',
                             cell_type_col = 'sctype_id') {
  Idents(seurat) <- cluster_col
  Labels <- table(seurat[[cell_type_col]], seurat@active.ident) %>% 
    as.data.frame() %>% 
    filter(Var1 != "Unknown") %>% 
    group_by(Var2) %>%
    filter(Freq == max(Freq)) %>%
    as.data.frame()
  
  colnames(Labels) <- c("Cell","seurat_clusters","Freq")
  
  return (Labels)
}
