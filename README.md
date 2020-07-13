# GNET2
Gene regulatory network modeling tool (GNET2)

[Links on Bioconductor](https://bioconductor.org/packages/devel/bioc/html/GNET2.html)

License
-------
© Contributors, 2019. Licensed under an [Apache-2] license.

Contribute to GNET2
---------------------
GNET2 has been developed and used by the Bioinformatics, Data Mining and Machine Learning Laboratory (BDM)
. Help from every community member is very valuable to make the package better for everyone.
Checkout the [Lab Page](http://calla.rnet.missouri.edu/cheng/).

## Install
(Recommended) To install our latest versions from GitHub, open R and type:
```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
library(devtools)
install_github("chrischen1/GNET2")
```

Older versions can be installed from Bioconductor:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GNET2", version = "3.9")
```

Vignette can be viewed with following code:

```
vignette('run_gnet2')
```

## Quick start
GNET2 is a R package used for build regulator network and cluster genes to functional groups with E-M process.
It iteratively perform TF assigning and Gene assigning, until the assignment of genes did not change, or max number of iterations is reached.

```
library(GNET2)
```

We first generate random expression data and a list of regulator gene names. 

The input is typically a p by n matrix of expression data of p genes and n samples, for example log2 RPKM from RNA-Seq.

```
set.seed(1)
gene_num = 1000
tf_num = 100
exp_data <- matrix(rnorm(gene_num*12),gene_num,12)
tf_list <- paste0('TF',1:tf_num)
rownames(exp_data) <- c(tf_list,paste0('gene',1:(nrow(exp_data)-length(tf_list))))
colnames(exp_data) <- paste0('sample_',1:ncol(exp_data))
se <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=exp_data))
```


The module construction process make take a few time, depending on the size of data and maximum iterations allowed.

```
gnet_result <- gnet(se,tf_list)
```

Show the total number of modules

```
print('Total number of modules:')
print(gnet_result$modules_count)
```

Show the regulatory genes and target genes in the first module

```
group_index = 1
print('Regulators in module 1:')
print(gnet_result$regulators[[group_index]])
print('Targets in module 1:')
print(gnet_result$target_genes[[group_index]])
```

Return the interations and their scores as adjacency matrix

```
mat <- extract_edges(gnet_result)
print(dim(mat))
```

Plot the tree of the first group

```
plot_gene_group(gnet_result,group_idx = group_index)
```

Plot the regulators module and heatmap of the expression inferred
downstream genes for each sample. It can be interpreted as two parts:
the bars at the top shows how samples are splited by the regression
tree and the heatmap at the bottom shows how downstream genes are
regulated by each subgroup determined by the regulators.

```
plot_tree(gnet_result,group_idx = group_index)
```


Plot the correlation of each group and auto detected knee point.

```
group_above_kn <- plot_group_correlation(gnet_result)
print(group_above_kn)
```

The group indices in group_above_kn can been seen as a list of indices of the data point with correlation higher than the knee point, and you may consider use them only for further analysis.
