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
CreateCRscript_Count <- function(file = NULL,
                                 cluster_base_directory,
                                 id_prefix,
                                 sample_prefix,
                                 sample_numbers,
                                 expect_cells,
                                 transcriptome){
  if (!(is.null(file))){
    #load file parameters here
  } else {
    #parse sample fastq numbers
  }
  if (str_detect(sample_numbers, ",")){
    first_num <- unlist(str_split(s, pattern = ","))[1]
    last_num <- unlist(str_split(s, pattern = ","))[-1]
    final_sample_numbers <- seq(first_num, last_num)
  } else {
    final_sample_numbers <- sample_numbers
  }
  for (n in sample_numbers){
    final_id = str_glue("{id_prefix}_{n}")
    filename = str_glue("{final_id}.pbs")
    
    write(
      x = str_glue(
        "#PBS -N {final_id} 
#PBS -l select=1:ncpus=16:mem=60gb:interconnect=1g,walltime=48:00:00
#PBS -m abe\n#PBS -M grangerb@musc.edu\n
source ~/.bashrc\n
#creating run directory
mkdir -p {cluster_base_directory}/run
cd {cluster_base_directory}/run\n
cellranger count --id={final_id} \\ 
--transcriptome=/scratch1/bryangranger/transcriptome/{transcriptome_reference} \\ 
--fastqs={cluster_base_directory}/fastq \\ 
--sample={sample_prefix}-{n} \\ 
--expect-cells={expect_cells} \n

#cellbender
conda activate cellbender\n
#creating cellbender output directory
mkdir -p {cluster_base_directory}/cellbender_out
cd {cluster_base_directory}/cellbender_out\n
cellbender remove-background \\ 
--input {cluster_base_directory}/run/{final_id}/outs/raw_feature_bc_matrix.h5 \\ 
--output {final_id}_cellbender_matrix.h5 \\ 
--expected-cells {expect_cells} \\ 
--total_droplets_included XXXXX \\ 
--fpr 0.01 \\
--epochs 150 
"),
      file = filename, 
      append = TRUE)
  }
}
  