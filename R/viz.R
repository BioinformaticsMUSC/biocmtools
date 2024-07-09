#' Create a proportional bar plot
#'
#' @param object The Seurat object
#' @param x_axis The column in the Seurat object meta.data to use as the x axis
#' @param split_by The column in the Seurat object meta.data used to split each bar
#' @return NULL
#' @importFrom dplyr %>% add_count select distinct mutate
#' @importFrom ggplot2 ggplot geom_bar ggtitle
#' @importFrom stringr str_detect
#' @examples
#' BarProp(seurat_obj, 'seurat_clusters', 'cell_type')
#' @export
BarProp <- function(seurat, x_axis, split_by, add_counts = FALSE) {
  
  df <- seurat@meta.data %>%
    as.data.frame() %>%
    add_count(!!as.name(x_axis), name = 'total_cluster_cells') %>%
    add_count(!!as.name(x_axis), !!as.name(split_by), name = 'split_counts') %>%
    select(!!as.name(x_axis), !!as.name(split_by), total_cluster_cells, split_counts) %>%
    distinct() %>%
    mutate(Percentage = split_counts / total_cluster_cells)
  
    fig <- ggplot(df, aes(x = factor(!!as.name(x_axis), levels = sort(unique(!!as.name(x_axis)))),
               y = Percentage,
               fill = !!as.name(split_by),
               col = !!as.name(split_by))) +
    geom_bar(stat = "identity") + ggtitle(str_glue("Bar Proportion Plot")) +
    theme(axis.text.x = element_text(angle = 60, vjust = 0.9, hjust = 1)) +
      xlab(x_axis)
    
    if (isTRUE(add_counts)) {
      fig <- fig +
        geom_text(aes(label = split_counts), vjust = 0.9, color = "black")
    }
    return (fig)
}

#' Create a list of colors for a palette
#'
#' @param n Number of colors
#' @return A named list of colors to use in plotting.
#' @importFrom scales hue_pal
#' @examples
#' get_hue_pal_list(12)
#' @export
get_hue_pal_list <- function (n) {
  pal <- scales::hue_pal()(n)
  color_list <- list()
  for (i in 1:n) {
    color_list[[paste0(i-1)]] <- pal[i]
  }
  return (color_list)
}
#' Create a dimplot with markers
#'
#' @param object The Seurat object
#' @param group_by The column in the Seurat object meta.data to group data by (e.g. "Cell" or "seurat_clusters")
#' @param split_by The column in the Seurat object meta.data used to split each bar
#' @param n The number of markers per identity to display
#' @return A plot object featuring a DimPlot next to a table of markers
#' @importFrom scales hue_pal
#' @importFrom presto wilcoxauc top_markers
#' @importFrom Seurat DimPlot
#' @import patchwork
#' @examples
#' MarkerDimPlot(seurat_obj, 'cell_type')
#' @export
MarkerDimPlot <- function (object, group_by, n = 20){
  #create markers table
  cat('Finding markers...\n')
  markers <- presto::wilcoxauc(object, group_by)
  top <- presto::top_markers(markers, n = n)
  if ("10" %in% colnames(top)) {
    rank_i <- which(colnames(top) == 'rank')
    sorted_colnames <- c("rank", sort(as.numeric(colnames(top[,-c(rank_i)]))))
  } else {
    sorted_colnames <- colnames(top)
  }
  
  top_df <- top[,sorted_colnames]
  
  cat('done!\nNow creating DimPlot...')
  p1 <- DimPlot(seurat, group.by = group_by)
  p2 <- gridExtra::tableGrob(top_df)
  
  if (ncol(top_df) < 9) {
    (p1 | p2) + plot_annotation(title = "DimPlot with markers")
  } else {
    (p1 / p2) + plot_annotation(title = "DimPlot with markers") 
  }
  
}
#' Create a clean dimplot
#'
#' @param object The Seurat object
#' @param group_by The column in the Seurat object meta.data to use to group cells (cluster, sample, etc.)
#' @importFrom dplyr %>% group_by summarize left_join mutate
#' @importFrom ggplot2 ggplot
#' @importFrom ggrastr geom_point_rast 
#' @importFrom ggrepel geom_text_repel
#' @examples
#' CleanDimPlot(seurat_obj, 'seurat_clusters')
#' @export
CleanDimPlot <- function (object, group_by) {
  umap <- as.data.frame(Embeddings(object, reduction = "umap"))
  meta <- as.data.frame(object@meta.data)
  
  df <- cbind(umap,meta)%>% 
    as.data.frame() 
  
  if (!(is.factor(df[[group_by]]))){
    df[[group_by]] <- factor(df[[group_by]], levels = unique(df[[group_by]]))
  }
  
  
  label <- data.frame(join_col=levels(df[[group_by]]),label=levels(df[[group_by]]))
  
  label_2 <- df %>% 
    group_by(!!as.name(group_by)) %>% 
    summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)) %>% 
    mutate(join_col = !!as.name(group_by)) %>%
    left_join(label) %>%
    as.data.frame() 
  
  col <- get_hue_pal_list(length(levels(df[[group_by]])))

  ggplot(df, aes(x=UMAP_1, y=UMAP_2)) +
    ggrastr::geom_point_rast(aes(colour = .data[[group_by]]), size=0.5) +
    ggrepel::geom_text_repel(data = label_2, aes(label = label),
                             color = "black",
                             segment.colour = "grey60",
                             box.padding = unit(0.25, "lines"),
                             point.padding = unit(0.5, "lines"),
                             nudge_x = .15,
                             nudge_y = 1,
                             size = 6) + 
    scale_colour_manual(values = col)+
    theme_void() +
    theme(legend.position="none")
}
#' Create a clean dimplot split by a specific category
#'
#' @param object The Seurat object
#' @param group_by The column in the Seurat object meta.data to use to group cells (cluster, sample, etc.)
#' @param split_by The column in the Seurat object to split the DimPlot by
#' @param ncol Number of columns in final figure (ignored for 3 or fewer plots)
#' @param split_order Order that split DimPlots will appear
#' @importFrom dplyr %>% group_by summarize left_join mutate
#' @importFrom ggplot2 ggplot
#' @importFrom ggrastr geom_point_rast 
#' @importFrom ggrepel geom_text_repel
#' @importFrom patchwork wrap_plots plot_layout
#' @importFrom stringr str_replace
#' @examples
#' MultiDimPlot(seurat_obj, 'seurat_clusters', 'condition', ncol = 4)
#' @export
MultiDimPlot <- function(seurat, group_by, split_by, ncol = 2, split_order = NULL) {
  umap <- as.data.frame(Embeddings(seurat, reduction = "umap"))
  meta <- as.data.frame(seurat@meta.data)
  
  df <- cbind(umap,meta)%>% 
    as.data.frame() 
  
  if (!(is.factor(df[[group_by]]))){
    df[[group_by]] <- factor(df[[group_by]], levels = unique(df[[group_by]]))
  }
  
  if (!(is.factor(df[[split_by]]))){
    if (is.null(split_order)) {
      df[[split_by]] <- factor(df[[split_by]], levels = unique(df[[split_by]]))
    } else {
      df[[split_by]] <- factor(df[[split_by]], levels = split_order)
    }
  }
  
  label <- data.frame(join_col=levels(df[[group_by]]),label=levels(df[[group_by]]))
  
  label_2 <- df %>% 
    group_by(!!as.name(group_by)) %>% 
    summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)) %>% 
    mutate(join_col = !!as.name(group_by)) %>%
    left_join(label) %>%
    as.data.frame() 
  
  col <- biocmtools::get_hue_pal_list(length(levels(df[[group_by]])))
  new_col = list()
  i = 0
  for (f in levels(df[[group_by]])) {
    new_col[[f]] = col[[as.character(i)]]
    i = 1 + i
  }
  viz=list()
  for (comp in levels(df[[split_by]])) {
    viz[[comp]] <- df %>%
      filter(!!as.name(split_by) == comp) %>%
      ggplot(aes(x=UMAP_1, y=UMAP_2,
                 color= .data[[group_by]])) +
      geom_point(size=0.5) + theme_void() +
      ggrepel::geom_text_repel(data = label_2, aes(label = label),
                               color = "black",
                               segment.colour = "grey60",
                               box.padding = unit(0.25, "lines"),
                               point.padding = unit(0.5, "lines"),
                               nudge_x = .15,
                               nudge_y = 1,
                               size = 4) +
      theme(legend.position = "none") +
      scale_color_manual(values = biocmtools::get_hue_pal_list(length(levels(df[[group_by]])))) +
      ggtitle(stringr::str_replace(comp, "_", " "))
    
  }
  
  if (length(levels(df[[split_by]])) < 4) {
    wrap_plots(viz) + plot_layout(nrow = 1)
  } else {
    wrap_plots(viz) + plot_layout(ncol = ncol)
  }
}

#' Create a viz with a UMAP and a bar plot showing percentage of cells assigned to clusters
#'
#' @param object The Seurat object
#' @param cluster_col The column in the Seurat object meta.data to use to determine clusters (e.g. seurat_clusters)
#' @param condition_col The column in the Seurat object denoting condition
#' @param replicate_col The column in the Seurat object denoting replicates (male, female, etc.) or cell types
#' @param data If TRUE, returns a dataframe
#' @importFrom dplyr %>% group_by select add_count mutate distinct rename
#' @importFrom tidyr pivot_wider unite
#' @import ggplot2
#' @importFrom ggrastr geom_point_rast 
#' @importFrom ggrepel geom_text_repel
#' @import patchwork 
#' @importFrom stringr str_glue
#' @importFrom Seurat DimPlot
#' @examples
#' ClusterProp(seurat_obj, cluster_col = 'seurat_clusters', condition_col = 'activation', replicate_col = 'tcell_type' = 4)
#' @export
ClusterProp <- function(seurat, cluster_col, condition_col, replicate_col, data = FALSE,
                        save_dir = NULL) {
  if (is.null(save_dir)){
    save_dir = getwd()
  } 
  
  df <- seurat@meta.data %>%
    dplyr::select(!!as.name(condition_col), !!as.name(cluster_col), !!as.name(replicate_col)) %>%
    add_count(!!as.name(condition_col), name = "total_cells_in_condition") %>%
    group_by(!!as.name(cluster_col), !!as.name(condition_col)) %>%
    mutate(n_cells_in_cluster = n()) %>%
    mutate(pct = n_cells_in_cluster / total_cells_in_condition) %>%
    distinct() %>%
    dplyr::rename(cluster = !!as.name(cluster_col))
  
  p1 <- DimPlot(seurat, group.by = cluster_col, label = T) + ggtitle(str_glue('UMAP - {cluster_col}')) +
    scale_fill_manual(values = biocmtools::get_hue_pal_list(length(unique(df$cluster))),
                      drop = TRUE) +
    scale_color_manual(values = biocmtools::get_hue_pal_list(length(unique(df$cluster))),
                       drop = TRUE)
  for (tcell in unique(df[[replicate_col]])) {
    t_df <- df %>% filter(!!as.name(replicate_col) == tcell)
    for (act in unique(t_df[[condition_col]])) {
      p2 <- ggplot(t_df %>% filter(!!as.name(condition_col) == act), 
                   aes(x = factor(cluster, levels = seq(0, length(unique(df$cluster))-1)), 
                       y = pct,
                       color = cluster,
                       fill = cluster)) + 
        geom_bar(stat = 'identity') +
        xlab('Cluster') +
        ylab(str_glue('Percentage of {act} cells by cluster')) +
        labs(fill = 'Cluster', color = "Cluster") +
        scale_fill_manual(values = biocmtools::get_hue_pal_list(length(unique(df$cluster))),
                          drop = TRUE) +
        scale_color_manual(values = biocmtools::get_hue_pal_list(length(unique(df$cluster))),
                           drop = TRUE)
      
      
      (p1 | p2) + plot_annotation(title = str_glue('Percentage of {act} {tcell} cells by cluster'))
      
      ggsave(file.path(save_dir, str_glue('{tcell}_{act}_prop_per_cluster.pdf')), width = 10, height = 8)
    }
  }
  if (isTRUE(data)) {
    new_df <- df %>%
      tidyr::unite("combos", !!(as.name(replicate_col)), !!as.name(condition_col), remove = FALSE) %>%
      tidyr::pivot_wider(names_from = combos,
                         id_cols = cluster,
                         values_from = c(pct, n_cells_in_cluster)) %>%
      arrange(factor(cluster, levels = seq(0, length(cluster)-1)))
    return (new_df)
      }
   
}
