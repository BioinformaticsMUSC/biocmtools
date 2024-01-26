#' Visualize the amount of cells removed due to pMito filtering
#'
#' @param seurat_obj The Seurat object
#' @param pMito_thresh The label (celltype, cluster, etc) column for organizing groups
#' @param group_col The condition column (treatment, etc) used for comparing between
#' @param pMito_col The background group for use in the tests (usually "control" or similar)
#' @return ggplot
#' @importFrom dplyr mutate group_by add_count distinct case_when 
#' @importFrom ggplot2 ggplot aes geom_bar coord_flip geom_text ylab scale_fill_manual labs
#' @importFrom stringr str_glue
#' @examples
#' FilterVis_pMito(seurat_obj, pMito_thresh = 15, group_col = "sample", pMito_col = "pMito")
#' @export
FilterVis_pMito <- function(seurat_obj, 
                            pMito_thresh, 
                            group_col, 
                            pMito_col = "pMito") {
  tmp <- seurat_obj@meta.data |>
    dplyr::mutate(pmito_filter = case_when(!!as.name(pMito_col) > pMito_thresh ~ "Remove", .default = "Keep"))
  
  tmp |> 
    dplyr::group_by(!!as.name(group_col)) |> 
    dplyr::add_count(pmito_filter) |> 
    dplyr::distinct(!!as.name(group_col), pmito_filter, .keep_all = T) |>
    dplyr::mutate(extra = ifelse(pmito_filter == 'Remove', -n, n)) |>
    ggplot2::ggplot(ggplot2::aes(
      x=sample,
      y=extra,
      fill=pmito_filter
    )) +
    ggplot2::geom_bar(stat='identity') +
    ggplot2::coord_flip() +
    ggplot2::geom_text(aes(label=n, hjust=0.5+sign(extra)/2+sign(extra)*0.2), color='white') +
    ggplot2::ylab("Cells") +
    ggplot2::scale_fill_manual(values = c("blue", "red")) +
    ggplot2::labs(fill = "pMito Filtering",
                  title = "Cells removed via pMito filtering",
                  subtitle = stringr::str_glue("Threshold: {pMito_thresh}%"))
}

# filter_vis_pmito(big_seurat, pMito_thresh=15, group_col = "sample")
# filter_vis_pmito(unfilt, pMito_thresh=5, group_col = "sample")

#' Visualize the amount of cells removed due to nUMI filtering
#'
#' @param seurat_obj The Seurat object
#' @param umi_thresh The label (celltype, cluster, etc) column for organizing groups
#' @param group_col The condition column (treatment, etc) used for comparing between
#' @param umi_col The background group for use in the tests (usually "control" or similar)
#' @return ggplot
#' @importFrom dplyr mutate group_by add_count distinct case_when 
#' @importFrom ggplot2 ggplot aes geom_bar coord_flip geom_text ylab scale_fill_manual labs
#' @importFrom stringr str_glue
#' @examples
#' FilterVis_nUMI(seurat_obj, umi_thresh = 20000, group_col = "sample", umi_col = "nUMI")
#' @export
FilterVis_nUMI <- function(seurat, umi_thresh, group_col, umi_col = "nUMI") {
  tmp <- seurat@meta.data |>
    mutate(umi_filter = case_when(!!as.name(umi_col) > umi_thresh ~ "Remove", .default = "Keep"))
  
  tmp$umi_filter <- factor(tmp$umi_filter, levels = c("Keep", "Remove"))
  
  tmp |> 
    group_by(!!as.name(group_col)) |> 
    add_count(umi_filter) |> 
    distinct(!!as.name(group_col), umi_filter, .keep_all = T) |>
    mutate(extra = ifelse(umi_filter == 'Remove', -n, n)) |>
    ggplot(aes(
      x=sample,
      y=extra,
      fill=umi_filter
    )) +
    geom_bar(stat='identity') +
    coord_flip() +
    geom_text(aes(label=n, hjust=0.5+sign(extra)/2+sign(extra)*0.2), color='white') +
    ylab("Cells") +
    scale_fill_manual(values = c("blue", "red")) +
    labs(fill = "nUMI Filtering",
         title = "Cells removed via UMI filtering",
         subtitle = stringr::str_glue("Threshold: {umi_thresh} UMI"))
}