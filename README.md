# GNET2
Gene regulatory network modeling tool (GNET2)

License
-------
© Contributors, 2019. Licensed under an [Apache-2] license.

Contribute to GNET2
---------------------
GNET2 has been developed and used by the Bioinformatics, Data Mining and Machine Learning Laboratory (BDM)
. Help from every community member is very valuable to make the package better for everyone.
Checkout the [Lab Page](http://calla.rnet.missouri.edu/cheng/).

## Install
To install, open R and type:
```
install.packages("devtools")
library("devtools")
install_github("chrischen1/GNET2")
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
rnaseq_data <- matrix(rnorm(gene_num*12),gene_num,12)
tf_list <- paste0('TF',1:tf_num)
rownames(rnaseq_data) <- c(tf_list,paste0('gene',1:(nrow(rnaseq_data)-length(tf_list))))
colnames(rnaseq_data) <- paste0('sample_',1:ncol(rnaseq_data))

head(rnaseq_data)
```


The module construction process make take a few time, depending on the size of data and maximum iterations allowed.

```
gnet_result = gnet(rnaseq_data,tf_list)
```

Plot the tree of the first group

```
plot_gene_group(gnet_result,group_idx = 11)
```

Plot the regulators module and heatmap of the expression inferred
downstream genes for each sample. It can be interpreted as two parts:
the bars at the top shows how samples are splited by the regression
tree and the heatmap at the bottom shows how downstream genes are
regulated by each subgroup determined by the regulators.

```
plot_tree(gnet_result,group_idx = 0)
```


Plot the correlation of each group and auto detected knee point.

```
group_above_kn <- plot_group_correlation(gnet_result)
print(group_above_kn)
```

The group indices in group_above_kn can been seen as a list of indices of the data point with correlation higher than the knee point, and you may consider use them only for further analysis.