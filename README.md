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
<ul>
	<li>seurat = seurat object</li>
	<li>x_axis = the X axis of the plot, usually cell type or Seurat clusters</li>
	<li>split_by = the category to split each bar by, usually sample, condition, or any other categorical metadata column</li>
</ul>


```
BarProp(seurat, x_axis = 'Cell', split_by = 'sample')
```

### MarkerDimPlot

This function creates a plot that includes a DimPlot from Seurat along with a table of markers. This is useful to quickly visualize markers per cluster or cell type. It requires the following parameters:
<ul>
	<li>seurat = seurat object</li>
	<li>group_by = the category to show in the DimPlot, usually 'Cell', 'Seurat Clusters', or another resolution</li>
	<li>n = the number of markers to show for each cluster or cell type. The default is 20.</li>
</ul>

This function uses presto to generate markers and patchwork to combine the plot with the table. The resulting plot may need a large saved image size (width = 18 inches, height = 18 inches for example).

```
MarkerDimPlot(seurat, group_by = 'Cell', n = 20)
```

### CleanDimPlot

This function creates a clean DimPlot from a Seurat object, without axes or a legend. It takes the following parameters:
- seurat = seurat object
- group_by = category in Seurat meta data to group the cells by

```
CleanDimPlot(seurat, group_by = 'Cell')
```

### MultiDimPlot

This function creates an array of Clean DimPlots split by a specified category in the Seurat meta data. It uses the following parameters:
- seurat = seurat object
- group_by = category in Seurat meta data to group the cells by
- split_by = category in Seurat meta data to split the plots by
- ncol = number of columns for the final figure (ignored for 3 or fewer plots)
- split_order = order that split plots will appear

```
MultiDimPlot(seurat, group_by = 'Cell', split_by = 'condition', ncol = 4)
```

