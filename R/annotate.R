# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/blob/master/LICENSE)
# Written by Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi>, June 2021

#' Prepare gene sets for sc-type
#'
#' @param path_to_db_file URL to excel sheet with gene sets
#' @param cell_type Type of cell for gene sets (e.g. "Immune system")
#' @return list
#' @import dplyr
#' @import openxlsx
#' @import HGNChelper
#' @examples
#' gene_sets_prepare("file_url", "Immune system")
#' @export
gene_sets_prepare <- function(path_to_db_file, cell_type){
  
  cell_markers = openxlsx::read.xlsx(path_to_db_file)
  cell_markers = cell_markers[cell_markers$tissueType == cell_type,] 
  cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1); cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  
  # correct gene symbols from the given DB (up-genes)
  cell_markers$geneSymbolmore1 = sapply(1:nrow(cell_markers), function(i){
    
    markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore1[i],",")))
    markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(markers_all)
    
    if(length(markers_all) > 0){
      suppressMessages({markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))})
      paste0(markers_all, collapse=",")
    } else {
      ""
    }
  })
  
  # correct gene symbols from the given DB (down-genes)
  cell_markers$geneSymbolmore2 = sapply(1:nrow(cell_markers), function(i){
    
    markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore2[i],",")))
    markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(markers_all)
    
    if(length(markers_all) > 0){
      suppressMessages({markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))})
      paste0(markers_all, collapse=",")
    } else {
      ""
    }
  })
  
  cell_markers$geneSymbolmore1 = gsub("///",",",cell_markers$geneSymbolmore1);cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 = gsub("///",",",cell_markers$geneSymbolmore2);cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  
  gs = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore1[j]),",")))); names(gs) = cell_markers$cellName
  gs2 = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore2[j]),",")))); names(gs2) = cell_markers$cellName
  
  return (list(gs_positive = gs, gs_negative = gs2))
}
# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/blob/master/LICENSE)
# Written by Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi>, June 2021
#
#' Calculate sc-type score
#'
#' @param scRNAseqData input scRNA-seq matrix (rownames - genes, column names - cells)
#' @param scaled indicates whether the matrix is scaled (TRUE by default)
#' @param gs list of gene sets positively expressed in the cell type
#' @param gs2 list of gene sets that should not be expressed in the cell type (NULL if not applicable)
#' @return data.frame
#' @importFrom Seurat DimPlot
#' @importFrom ggplot2 ggtitle
#' @import dplyr
#' @import HGNChelper
#' @examples
#' sctype_score("file_url", "Immune system")
#' @export
sctype_score <- function(scRNAseqData, scaled = !0, gs, gs2 = NULL, gene_names_to_uppercase = !0, ...){
  
  # check input matrix
  if(!is.matrix(scRNAseqData)){
    warning("scRNAseqData doesn't seem to be a matrix")
  } else {
    if(sum(dim(scRNAseqData))==0){
      warning("The dimension of input scRNAseqData matrix equals to 0, is it an empty matrix?")
    }
  }
  
  # marker sensitivity
  marker_stat = sort(table(unlist(gs)), decreasing = T); 
  marker_sensitivity = data.frame(score_marker_sensitivity = scales::rescale(as.numeric(marker_stat), to = c(0,1), from = c(length(gs),1)),
                                  gene_ = names(marker_stat), stringsAsFactors = !1)
  
  # convert gene names to Uppercase
  if(gene_names_to_uppercase){
    rownames(scRNAseqData) = toupper(rownames(scRNAseqData));
  }
  
  # subselect genes only found in data
  names_gs_cp = names(gs); names_gs_2_cp = names(gs2);
  gs = lapply(1:length(gs), function(d_){ 
    GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
  gs2 = lapply(1:length(gs2), function(d_){ 
    GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs2[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
  names(gs) = names_gs_cp; names(gs2) = names_gs_2_cp;
  cell_markers_genes_score = marker_sensitivity[marker_sensitivity$gene_ %in% unique(unlist(gs)),]
  
  # z-scale if not
  if(!scaled) Z <- t(scale(t(scRNAseqData))) else Z <- scRNAseqData
  
  # multiple by marker sensitivity
  for(jj in 1:nrow(cell_markers_genes_score)){
    Z[cell_markers_genes_score[jj,"gene_"], ] = Z[cell_markers_genes_score[jj,"gene_"], ] * cell_markers_genes_score[jj, "score_marker_sensitivity"]
  }
  
  # subselect only with marker genes
  Z = Z[unique(c(unlist(gs),unlist(gs2))), ]
  
  # combine scores
  es = do.call("rbind", lapply(names(gs), function(gss_){ 
    sapply(1:ncol(Z), function(j) {
      gs_z = Z[gs[[gss_]], j]; gz_2 = Z[gs2[[gss_]], j] * -1
      sum_t1 = (sum(gs_z) / sqrt(length(gs_z))); sum_t2 = sum(gz_2) / sqrt(length(gz_2));
      if(is.na(sum_t2)){
        sum_t2 = 0;
      }
      sum_t1 + sum_t2
    })
  })) 
  
  dimnames(es) = list(names(gs), colnames(Z))
  es.max <- es[!apply(is.na(es) | es == "", 1, all),] # remove na rows
  
  return (es.max)
}

#' Assess batch correction by Harmony by printing out side by side dim plots before and after correction
#'
#' @param seurat The Seurat object
#' @param assay The assay to use, should be 'RNA' or 'SCT'
#' @param cell_type The column in the Seurat object that distinguishes replicates (eg. "patient")
#' @return Seurat object
#' @import Seurat 
#' @import plyr 
#' @examples
#' RunSCtype(seurat_obj, "Immune system")
#' @export
RunSCtype <- function(seurat, cell_type){
  
  cat('Running sc-type on seurat_clusters:\n')
  cat(levels(seurat$seurat_clusters))
  gs_list = gene_sets_prepare("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx", cell_type)
  
  es.max = sctype_score(scRNAseqData = seurat[[assay]]@scale.data, scaled = TRUE,
                        gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
  
  cL_results = do.call("rbind", lapply(unique(seurat@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(seurat@meta.data[seurat@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat@meta.data$seurat_clusters==cl)), 10)
  }))
  sctype_scores = cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
  
  if (length(unique(sctype$cluster)) != unique(sctype$cluster)) {
    warning("\nYou may have duplicate cluster labels - run sctypes manually to address this.")
  }
  
  sctype_scores <- distinct(sctype_scores, cluster, .keep_all = TRUE)
  
  cat("\nAdding cell types to Seurat object...")
  Idents(seurat) <- 'seurat_clusters'
  
  seurat@meta.data$sctype_id <- seurat@active.ident
  cat("done!\n")
  
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
