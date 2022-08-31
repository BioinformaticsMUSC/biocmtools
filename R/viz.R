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
BarProp <- function(seurat, x_axis, split_by) {
  
  df <- seurat@meta.data %>%
    as.data.frame() %>%
    add_count(!!as.name(x_axis), name = 'total_cluster_cells') %>%
    add_count(!!as.name(x_axis), !!as.name(split_by), name = 'split_counts') %>%
    select(!!as.name(x_axis), !!as.name(split_by), total_cluster_cells, split_counts) %>%
    distinct() %>%
    mutate(Percentage = split_counts / total_cluster_cells)
  
 #if (all(str_detect(df %>% select(!!as.name(x_axis)), 
  #                   pattern = "[:digit:]"))) {
  #  df <- df %>%
  #    mutate(!!as.name(x_axis) = as.numeric(!!as.name(x_axis)))
  #}
  
    ggplot(df, aes(x = !!as.name(x_axis),
               y = Percentage,
               fill = !!as.name(split_by),
               col = !!as.name(split_by))) +
    geom_bar(stat = "identity") + ggtitle(str_glue("Bar Proportion Plot")) +
    theme(axis.text.x = element_text(angle = 60, vjust = 0.9, hjust = 1))
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
