BiocManager::install('recount')

## Browse the vignetets for a quick description of how to use the package
library('recount')
project_info <- abstract_search('GSE32465')

## Download the gene-level RangedSummarizedExperiment data
download_study(project_info$project)

## Load the data
load(file.path(project_info$project, 'rse_gene.Rdata'))

## Browse the project at SRA
browse_study(project_info$project)

## View GEO ids
colData(rse_gene)$geo_accession

## Extract the sample characteristics
geochar <- lapply(split(colData(rse_gene), seq_len(nrow(colData(rse_gene)))), geo_characteristics)

## Note that the information for this study is a little inconsistent, so we
## have to fix it.
geochar <- do.call(rbind, lapply(geochar, function(x) {
  if('cells' %in% colnames(x)) {
    colnames(x)[colnames(x) == 'cells'] <- 'cell.line'
    return(x)
  } else {
    return(x)
  }
}))

## We can now define some sample information to use
sample_info <- data.frame(
  run = colData(rse_gene)$run,
  group = ifelse(grepl('uninduced', colData(rse_gene)$title), 'uninduced', 'induced'),
  gene_target = sapply(colData(rse_gene)$title, function(x) { strsplit(strsplit(x,
                                                                                'targeting ')[[1]][2], ' gene')[[1]][1] }),
  cell.line = geochar$cell.line
)

## Scale counts by taking into account the total coverage per sample
rse <- scale_counts(rse_gene)

## Add sample information for DE analysis
colData(rse)$group <- sample_info$group
colData(rse)$gene_target <- sample_info$gene_target


