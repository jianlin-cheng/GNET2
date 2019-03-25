set.seed(2019)

# 20 TF, 100 genes,10 samples
rnaseq_data <- matrix(rnorm(120*10),120,10)

rownames(rnaseq_data) <- c(paste0('tf',1:20),paste0('gene',1:100))
write.csv(rnaseq_data,paste0('~/Dropbox/MU/workspace/new/data/test/rnaseq_data.csv'))
write.table(paste0('tf',1:20),paste0('~/Dropbox/MU/workspace/new/data/test/tf_list.csv'),col.names = F,row.names = F,quote = F)

set.seed(2019)

rnaseq_data = matrix(rnorm(7000*20),7000,20)
colnames(rnaseq_data) <- paste0('sample_',1:ncol(rnaseq_data))
tf_list = paste0('tf_',1:500)
rownames(rnaseq_data) <- c(tf_list,paste0('gene_',1:6500))

write.csv(rnaseq_data,'~/Dropbox/MU/workspace/new/data/test2/rnaseq_data.csv',quote = F)
write.table(tf_list,'~/Dropbox/MU/workspace/new/data/test2/tf_list.csv',col.names = F,row.names = F,quote = F)