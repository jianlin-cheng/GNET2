## ----setup,echo = FALSE, include = TRUE----------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load_packages, message=FALSE, warning=FALSE, include=FALSE----------
library(GNET2)

## ----echo = TRUE---------------------------------------------------------
set.seed(1)
gene_num = 1000
tf_num = 100
rnaseq_data <- matrix(rnorm(gene_num*12),gene_num,12)
tf_list <- paste0('TF',1:tf_num)
rownames(rnaseq_data) <- c(tf_list,paste0('gene',1:(nrow(rnaseq_data)-length(tf_list))))
colnames(rnaseq_data) <- paste0('sample_',1:ncol(rnaseq_data))

head(rnaseq_data)

## ---- echo=TRUE----------------------------------------------------------
gnet_result = gnet(rnaseq_data,tf_list)

## ---- echo=TRUE,fig.width=10, fig.height=12------------------------------
plot_gene_group(gnet_result,group_idx = 11)

## ---- echo=TRUE,fig.width=8, fig.height=8--------------------------------
plot_tree(gnet_result,group_idx = 0)

## ---- echo=TRUE,fig.width=8, fig.height=8--------------------------------
group_above_kn <- plot_group_correlation(gnet_result)
print(group_above_kn)

