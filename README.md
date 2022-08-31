# biocmtools
## A collection of tools for the BioCM

This is a package for tools used by the BioCM.

## Installation

Biocmtools currently requires the following packages:

```
Seurat
ggplot2
patchwork
dplyr
stringr
scales
presto
```

To install biocmtools, please do so via devtools:


```
devtools::install_github("BioinformaticsMUSC/biocmtools")
```

## Usage

There are currently two main tools.

### BarProp

This function seamlessly creates a bar proportion plot from a Seurat object. It's parameters are:
seurat = seurat object
x_axis = the X axis of the plot, usually cell type or Seurat clusters
split_by = the category to split each bar by, usually sample, condition, or any other categorical metadata column

```
BarProp(seurat, x_axis = 'Cell', split_by = 'sample')
```

### MarkerDimPlot

This function creates a plot that includes a DimPlot from Seurat along with a table of markers. This is useful to quickly visualize markers per cluster or cell type. It requires the following parameters:
seurat = seurat object
group_by = the category to show in the DimPlot, usually 'Cell', 'Seurat Clusters', or another resolution
n = the number of markers to show for each cluster or cell type. The default is 20.

This function uses presto to generate markers and patchwork to combine the plot with the table. The resulting plot may need a large saved image size (width = 18 inches, height = 18 inches for example).

```
MarkerDimPlot(seurat, group_by = 'Cell', n = 20)
```

