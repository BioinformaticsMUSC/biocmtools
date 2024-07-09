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

#' Quickly create volcano plots for each celltype
#'
#' @param de_results A table of differential expression results grouped by celltype
#' @param celltype_col The label column for organizing groups 
#' @param gene_col The column of genes
#' @param fold_change_col The name of the log fold change column
#' @param save Save individual plots for each celltype
#' @return List of ggplots
#' @importFrom dplyr mutate case_when filter 
#' @importFrom tidyr unite
#' @importFrom ggplot2 xlab ylab geom_vline geom_hline theme arrow annotate xlim ggsave ggtitle
#' @importFrom ggpubr ggscatter
#' @importFrom ggrepel geom_text_repel
#' @importFrom stringr str_glue
#' @examples
#' QuickVolcano(de_results, celltype_col = "celltype", gene_col = "gene", fold_change_col = "avg_log2FC", save = TRUE)
#' @export
QuickVolcano <- function(de_results, 
                         celltype_col, 
                         gene_col, 
                         fold_change_col,
                         p_val_adj,
                         test_group_name = NULL, save = TRUE) {
  
  if (!is.null(test_group_name)) {
    upreg_msg = paste0("Upreg in ", test_group_name)
    downreg_msg = paste0("Downreg in ", test_group_name)
  } else {
    upreg_msg = "Upregulated"
    downreg_msg = "Downregulated"
  }

  data <- de_results |>
    dplyr::mutate(LOG = -log10(!!as.name(p_val_adj))) |>
    dplyr::mutate(direction = dplyr::case_when(!!as.name(fold_change_col) > 0 ~ "Upreg",
                                               .default = "Downreg")) |>
    dplyr::mutate(significance = dplyr::case_when(abs(!!as.name(fold_change_col)) > 0.3 & !!as.name(p_val_adj) < 0.5 ~ "sig", 
                                                  .default = 'not_sig')) |>
    tidyr::unite("direction_sig", direction, significance, remove = FALSE)
  
  data$direction_sig <- factor(data$direction_sig, levels = c("Upreg_sig", "Upreg_not_sig", "Downreg_not_sig", "Downreg_sig"))

  max_x = round(max(abs(data[[fold_change_col]])), 0)
  
  
  plot_list = list()
  for (ct in unique(data[[celltype_col]])) {
    tmp_data <- data |> dplyr::filter(!!as.name(celltype_col) == ct)
    downreg_count = sum(tmp_data[[fold_change_col]] < 0)
    upreg_count = sum(tmp_data[[fold_change_col]] > 0)
    max_y <- max(dplyr::filter(tmp_data, !!as.name(p_val_adj) != 0)$LOG)
    plot_list[[paste0("c",ct)]] <- ggpubr::ggscatter(tmp_data, 
                                 x = fold_change_col, 
                                 y = "LOG",
                                 color = "direction_sig", 
                                 palette = list(Upreg_sig="red",
                                                Upreg_not_sig="gray", 
                                                Downreg_not_sig="gray", 
                                                Downreg_sig="blue"),
                                 size = 1,
                                 alpha=0.3,
                                 shape=19)+
      ggplot2::xlab("log2(Fold Change)")+ 
      ggplot2::ylab("-log10(FDR)")+
      ggplot2::geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) + 
      ggplot2::geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) + 
      ggplot2::geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) + 
      ggplot2::geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      ggplot2::annotate("segment", x = 1, xend = 2, y = max_y + 15, yend = max_y + 15,
               size = 0.4, arrow = ggplot2::arrow(length = unit(0.1, "inches"), type = "closed")) + 
      ggplot2::annotate("segment", x = -1, xend = -2, y = max_y + 15, yend = max_y + 15,
               size = 0.4, arrow = ggplot2::arrow(length = unit(0.1, "inches"), type = "closed")) +
      ggplot2::annotate("text", x = .8 * max_x, y = max_y + 20, label = upreg_msg, size = 3) +
      ggplot2::annotate("text", x = .8 * -max_x, y = max_y + 20, label = downreg_msg, size = 3) +
      ggrepel::geom_text_repel(data = tmp_data, 
                      mapping = aes(label = !!as.name(gene_col)), 
                      size = 3,
                      box.padding = unit(0.4, "lines"),
                      point.padding = unit(0.4, "lines"))+
      #ggplot2::theme(legend.position="none") +
      ggplot2::xlim(-max_x,max_x) +
      ggplot2::ggtitle(ct) +
      ggplot2::labs(caption = stringr::str_glue("Downreg: {downreg_count}\nUpreg: {upreg_count}")) +
      ggplot2::theme(plot.caption = element_text(size=7),
                     legend.position="none")
    if (isTRUE(save)) {
      ggplot2::ggsave(stringr::str_glue("volc_{ct}.pdf"), width = 4, height = 4) 
    }
  }
  return (plot_list)
}
