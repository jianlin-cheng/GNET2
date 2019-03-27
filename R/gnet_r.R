#' Calculate Gaussian Likelihood score.
#' 
#' @param x A n by p matrix.
#' @param labels A vector of length n, indicating the group of rows.
#' @return The sum of log likelihood score of each group on each column.
#' @examples
#' calc_likelihood_score(x = matrix(rnorm(50*100),50,100), labels = c(rep(1,20),rep(2,30)))
#' @export
calc_likelihood_score = function(x,labels){
  score_total = 0
  for(i in unique(labels)){
    if(sum(labels==i) >1){
      x_i = x[labels==i,]
      for(j in 1:ncol(x_i)){
        x_i_j = x_i[,j]
        score_total = score_total+sum(dnorm(x_i_j,mean = mean(x_i_j),sd = sd(x_i_j),log = T))
      }
    }
  }
  return(score_total)
}

calc_likelihood_vec = function(x,labels){
  score_total = 0
  for(i in unique(labels)){
    if(sum(labels==i) >1){
      x_i = x[labels==i]
      score_total = score_total+sum(dnorm(x_i,mean = mean(x_i),sd = sd(x_i),log = T))
    }
  }
  return(score_total)
}

calc_correlation = function(x){
  if(nrow(x)<2 | ncol(x)<2){
    return(0)
  }else{
    x_cor = cor(x)
    return(mean(abs(x_cor[upper.tri(x_cor)])))
  }
}

#' Calculate Pearson correlation coefficient within each group.
#' 
#' @param x A n by p matrix.
#' @param labels A vector of length n, indicating the group of rows.
#' @return An array of Pearson correlation coefficient for each row, rows belong to the same group have same values.
#' @examples
#' get_correlation_list(x = matrix(rnorm(50*100),50,100), labels = c(rep(1,20),rep(2,30)))
#' @export
get_correlation_list = function(x,labels){
  cor_list = rep(0,length(labels))
  for(i in unique(labels)){
    if(sum(labels==i)==1){
      cor_list[labels==i] <- 0
    }else{
      cor_list[labels==i] = calc_correlation(x[,labels==i])
    }
  }
  return(cor_list)
}

get_leaf_labels = function(group_table,format_plot=F){
  if(!format_plot){
    group_table = group_table[,2:ncol(group_table),drop=F]
  }
  leaf_label = rep(0,ncol(group_table))
  next_label = 0
  for(i in 1:nrow(group_table)){
    current_label = group_table[i,]
    for(j in 0:1){
      leaf_label[current_label==j] = next_label
      next_label = next_label+1
    }
  }
  return(leaf_label)
}

#' R version of build regression tree based on Gaussian Likelihood score. Runs slower than the build_module() which is C++ version but gives flexibility of customization of R functions.
#' 
#' @param max_partition_level Maximum partition level in the tree.
#' @param cor_cutoff Cutoff for within group Pearson correlation coefficient, if all data belong to a node have average correlation greater or equal to this, the node would not split anymore.
#' @param min_divide_size Minimum number of data belong to a node allowed for further split of the node.
#' 
#' @return A matrix for sample informatrion for each partition level. First column is feature index used by the node and second is the value used to split, the rest of the columns are the split of sample: 0 means less or equal, 1 means greater and -1 means the sample does not belong to this node.
#' @examples
#' build_moduleR(X = matrix(rnorm(50*100),50,100), Y = matrix(rnorm(50*200),50,200))
#' @export
build_moduleR = function(X,Y,max_partition_level,cor_cutoff,min_divide_size){
  feature_remaining = 1:ncol(X)
  feature_num = 0
  subgroup_indicator = rep(0,nrow(X))
  groups = c(0)
  group_table = matrix(0, nrow = 0, ncol = nrow(X)+1)
  while (feature_num < max_partition_level & length(groups)>0 & length(feature_remaining)>0){
    best_score = -(10^6)
    best_feature = -1
    for(i in groups){
      current_split_group = subgroup_indicator == i
      y_current_split = Y[current_split_group,]
      for (j in feature_remaining){
        feature_vals = X[current_split_group,j]
        divide_vals = feature_vals[feature_vals != max(feature_vals)]
        score_nosplit = calc_likelihood_score(y_current_split,rep(1,length(feature_vals)))
        for (k in divide_vals){
          subgroup_divide = 1-(feature_vals <= k)
          score = calc_likelihood_score(y_current_split,subgroup_divide) - score_nosplit
          if (score > best_score){
            best_score = score
            best_feature = j
            subgroup_labels = rep(-1,length(current_split_group))
            subgroup_labels[current_split_group] = subgroup_divide
          }
        }
      }
    }
    group_table = rbind(group_table,c(best_feature-1,subgroup_labels))
    subgroup_indicator_new = get_leaf_labels(group_table)
    subgroup_indicator_new[subgroup_indicator==-1] = -1
    subgroup_indicator = subgroup_indicator_new
    for (i in 0:1){
      new_divide_i = subgroup_labels==i
      if (sum(new_divide_i) <= min_divide_size){
        subgroup_indicator[new_divide_i] = -1
      }else if(calc_correlation(t(Y[new_divide_i,])) >= cor_cutoff){
        subgroup_indicator[new_divide_i] = -1
      }
    }
    groups = unique(subgroup_indicator[subgroup_indicator!=-1])
    feature_remaining = feature_remaining[feature_remaining!=best_feature]
    feature_num = feature_num+1
  }
  return(group_table)
}

assign_tf = function(regulator_data,gene_data,gene_group_table,reg_names,
                             min_group_size,max_partition_level,cor_cutoff,min_divide_size){
  X = t(regulator_data)
  group_labels = unique(gene_group_table)
  group_labels = group_labels[group_labels!=-1]
  reg_group_table = matrix(0,nrow = 0,ncol = ncol(regulator_data))
  group_table = matrix(0,nrow = 0,ncol = ncol(regulator_data)+2)
  i = 0
  for (group_idx in group_labels){
    gene_idx = gene_group_table == group_idx
    if (sum(gene_idx) >= min_group_size){
      Y = t(gene_data[gene_idx,])
      group_table_i = build_module(X,Y,max_partition_level,cor_cutoff,min_divide_size)
      group_table_i = group_table_i[apply(group_table_i, 1, function(x)(sum(x!=0)))!=0,]
      reg_group_table = rbind(reg_group_table,get_leaf_labels(group_table_i))
      group_table = rbind(group_table,cbind(i,group_table_i))
      i = i+1
    }
  }
  return(list(group_table,reg_group_table))
}

assign_tfR = function(regulator_data,gene_data,gene_group_table,reg_names,
                             min_group_size,max_partition_level,cor_cutoff,min_divide_size){
  X = t(regulator_data)
  group_labels = 0:max(gene_group_table)
  reg_group_table = matrix(0,nrow = 0,ncol = ncol(regulator_data))
  group_table = matrix(0,nrow = 0,ncol = ncol(regulator_data)+2)
  i = 0
  for (group_idx in group_labels){
    gene_idx = gene_group_table == group_idx
    if (sum(gene_idx) >= min_group_size){
      Y = t(gene_data[gene_idx,])
      group_table_i = build_moduleR(X,Y,max_partition_level,cor_cutoff,min_divide_size)
      reg_group_table = rbind(reg_group_table,get_leaf_labels(group_table_i))
      group_table = rbind(group_table,cbind(i,group_table_i))
      i = i+1
    }
  }
  return(list(group_table,reg_group_table))
}

assign_gene = function(gene_data,reg_group_table){
  gene_group_table = rep(-1,nrow(gene_data))
  for(gene_idx in 1:nrow(gene_data)){
    exp_gene = as.numeric(gene_data[gene_idx,])
    gene_group_table[gene_idx] = which.max(
      apply(reg_group_table, 1, function(x)calc_likelihood_vec(exp_gene,x)))-1
  }
  return(gene_group_table)
}

#' Detect the knee point of the array.
#' 
#' @param vect A list of sorted numbers.
#' 
#' @return The index of the data point which is the knee.
#' @examples
#' kneepointDetection(sort(rnorm(100),T))
#' @export
kneepointDetection <-function (vect) {
  n <- length(vect)
  Vect <- vect
  a <- as.data.frame(cbind(1:n, Vect[1:n]))
  l <- lm(a[, 2] ~ a[, 1], data = a)
  MinError = 1e+08
  MinIndex = 1
  for (i in 2:(n - 2)) {
    a <- as.data.frame(cbind(1:i, Vect[1:i]))
    l1 <- lm(a[, 2] ~ a[, 1], data = a)
    e1 <- sum(abs(1 - a[, 2]))
    a <- as.data.frame(cbind((i + 1):n, Vect[(i + 1):n]))
    l <- lm(a[, 2] ~ a[, 1], data = a)
    l2 <- l
    e2 <- sum(abs(l$residuals))
    Error = e1 + e2
    if (MinError > Error) {
      MinError = Error
    }
      MinIndex = i
  }
  return(MinIndex)
}

run_gnet = function(gene_data,regulator_data,init_group_num = 5,max_partition_level = 3,
                    cor_cutoff = 0.9,min_divide_size = 3,min_group_size = 2,max_iter = 5){
  min_group_num=3
  message('Determine initial group number')
  reg_names = rownames(regulator_data)
  avg_cor_list <- c()
  for(i in 2:min(init_group_num,nrow(gene_data)-1)){
    gene_group_table = as.numeric(kmeans(gene_data,centers = i)$cluster)-1
    avg_cor_list <- c(avg_cor_list,mean(get_correlation_list(t(gene_data),gene_group_table)))
    if(length(avg_cor_list)>=min_group_num){
      o = order(avg_cor_list,decreasing = T)
      y = avg_cor_list[o]
      kn = kneepointDetection(y)
      groups_keep = o[1:max(kn,min_group_num)]
    }else{
      groups_keep = 1:length(avg_cor_list)
    }
    gene_group_table[!(gene_group_table %in% groups_keep)] <- -1
  }


  message('Building module networks')
  for (i in 1:max_iter) {
    message(paste('iteration',i))
    assign_reg_names = assign_tf(regulator_data,gene_data,gene_group_table,reg_names,
                                       min_group_size,max_partition_level,cor_cutoff,min_divide_size)
    gene_group_table_new = assign_gene(gene_data,assign_reg_names[[2]])
    if(all(length(gene_group_table)==length(gene_group_table_new)) &&
       all(gene_group_table==gene_group_table_new)){
      break
    }else{
      gene_group_table =  gene_group_table_new
    }
  }
  message('Done.')
  tree_table_all = assign_reg_names[[1]]
  colnames(tree_table_all) <- c('group','feature',colnames(gene_data))

  gene_list_all = NULL
  groups_unique = unique(gene_group_table)
  for (i in 1:length(groups_unique)) {
    gene_list_all = rbind(gene_list_all,
                           cbind.data.frame('gene'=rownames(gene_data)[gene_group_table==groups_unique[i]],
                                            'group'=i-1,stringsAsFactors=F))
  }
  return(list('reg_group_table'=tree_table_all,'gene_group_table'=gene_list_all))
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  # From http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    grid::grid.newpage()
    grid::pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#' Plot the regression tree given the index of a module.
#' 
#' @param gnet_result Results returned by gnet().
#' @param group_idx Index of the module..
#' 
#' @export
plot_tree <- function(gnet_result,group_idx){
  regulator_data = gnet_result$regulator_data
  gene_data = gnet_result$gene_data
  reg_group_table = gnet_result$reg_group_table
  gene_group_table = gnet_result$gene_group_table
  tree_i = reg_group_table[reg_group_table[,1]==group_idx,]
  label_list = c()
  for(i in 1:nrow(tree_i)){
    coverage_i = sum(tree_i[i,3:ncol(tree_i)]>=0)
    feature_name_i = rownames(regulator_data)[tree_i[i,2]+1]
    split_i = round(max(regulator_data[feature_name_i,tree_i[i,3:ncol(tree_i)]==0]),4)
    label_list = c(label_list,paste0(feature_name_i,'\n <= ',split_i,'\nCoverage: ',coverage_i))
  }
  leaf_idx = 0
  from_list <- to_list <- edge_labels <- c()
  for(i in 1:(nrow(tree_i)-1)){
    for(j in 0:1){
      found_node = F
      split_i_j = tree_i[i,3:ncol(tree_i)]==j
      for (k in (i+1):nrow(tree_i)) {
        node_k = tree_i[k,3:ncol(tree_i)]>=0
        if(identical(split_i_j,node_k)){
          from_list <- c(from_list,i)
          to_list <- c(to_list,k)
          edge_labels <- c(edge_labels,as.character(j==0))
          found_node = T
        }
      }
      if(!found_node){
        label_list <- c(label_list,paste('leaf',leaf_idx))
        from_list <- c(from_list,i)
        to_list <- c(to_list,length(label_list))
        edge_labels <- c(edge_labels,as.character(j==0))
        leaf_idx = leaf_idx+1
      }
    }
  }
  label_list <- c(label_list,paste('leaf',leaf_idx),paste('leaf',leaf_idx+1))
  from_list <- c(from_list,nrow(tree_i),nrow(tree_i))
  to_list <- c(to_list,length(label_list)-1,length(label_list))
  edge_labels <- c(edge_labels,'true','false')
  node_df <- create_node_df(n = length(label_list),type = "a",label = label_list,style = "filled",
                            color = "aqua",shape = "ellipse",width = 0.9,fontsize = 6)
  edge_df <-create_edge_df(from = from_list,to = to_list,label = tolower(edge_labels),fontsize = 6)
  graph <-create_graph(nodes_df = node_df,edges_df = edge_df,
                       attr_theme = NULL) %>% DiagrammeR::add_global_graph_attrs(attr_type = "graph", 
                      attr = c("layout", "rankdir"),
                      value = c("dot", "LR")) %>% DiagrammeR::add_global_graph_attrs(attr_type = "node", 
                      attr = c("fillcolor", "style", "fontname"), value = c("Azure","filled", "Helvetica"))
  render_graph(graph)
}

#' Plot the regulators module and heatmap of the expression inferred downstream genes for each sample. It can be interpreted as two parts: the bars at the top shows how samples are splited by the regression tree and the heatmap at the bottom shows how downstream genes are regulated by each subgroup determined by the regulators.
#' 
#' @param gnet_result Results returned by gnet().
#' @param group_idx Index of the module..
#' 
#' @export
plot_gene_group <- function(gnet_result,group_idx){
  gene_data <- gnet_result$gene_data
  regulator_data <- gnet_result$regulator_data
  reg_group_table <- gnet_result$reg_group_table
  gene_group_table <- gnet_result$gene_group_table

  exp_data1 = gene_data[gene_group_table$gene[gene_group_table$group==group_idx],]
  regulator_data1 = regulator_data[reg_group_table[reg_group_table[,1]==group_idx,2]+1,]
  group_table1 = reg_group_table[reg_group_table[,1]==group_idx,3:ncol(reg_group_table)]
  leaf_labels <- get_leaf_labels(group_table1,format_plot = T)
  row_order <- order(leaf_labels)
  group_table2 <- group_table1[,row_order]
  regulator_data2 <- regulator_data1[,row_order]
  exp_data2 <- exp_data1[,row_order]
  test_regulators_names <- rownames(regulator_data2)
  layout=matrix(c(1:length(test_regulators_names),rep(length(test_regulators_names)+1,
                                                      length(test_regulators_names)*2)),ncol=1)
  regulators_plist <- list()
  scaleFUN <- function(x) sprintf("%.3f", x)

  # add TF bars
  for(i in 1:length(test_regulators_names)){
    reg_data_mask <- group_table2[i,]==-1
    exp_val <- as.numeric(regulator_data2[i,])
    exp_val[reg_data_mask] <- NA
    lengend_low <- min(exp_val,na.rm = T)
    lengend_high <- max(exp_val,na.rm = T)
    exp_val1 <- rbind.data.frame(matrix(NA,nrow = 1,ncol = length(exp_val)),exp_val,stringsAsFactors=F)

    rownames(exp_val1) <- 1:nrow(exp_val1)
    exp_val.m <- reshape2::melt(exp_val1,id.vars = NULL)
    exp_val.m <- cbind.data.frame('y_idx'=rep(1:nrow(exp_val1),ncol(exp_val1)),exp_val.m,stringsAsFactors=F)
    exp_label = rep('',ncol(exp_val1))
    exp_label[group_table2[i,]==0] <- 'Low'
    exp_label[group_table2[i,]==1] <- 'High'

    p <- ggplot2::ggplot(exp_val.m, aes(variable, y_idx)) + geom_tile(aes(fill = value), colour = "white") +
      scale_x_discrete(labels=exp_label)+
      scale_fill_gradient(low = "darkgreen",high = "red",na.value = "white",
                          limits=c(lengend_low, lengend_high),
                          breaks=seq(lengend_low,lengend_high,length.out = 4),labels=scaleFUN)+
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.background = element_blank(),
            legend.title=element_blank(),panel.grid.minor = element_blank(),
            legend.key.size = unit(0.2, "cm"),
            axis.line = element_line(colour = "white"),legend.position="right",
            legend.box = "vertical",axis.title.x=element_blank(),
            axis.ticks.x=element_blank(),legend.text=element_text(size=7),
            axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
      labs(title = test_regulators_names[i])
    regulators_plist[[i]] <- p
  }
  # add heatmap
  exp_lengend_low <- min(exp_data2)
  exp_lengend_high <- max(exp_data2)
  test_data.m <- reshape2::melt(cbind.data.frame('gene'=rownames(exp_data2),exp_data2,stringsAsFactors=F),id.vars = 'gene')
  p <- ggplot2::ggplot(test_data.m, aes(variable, gene)) + geom_tile(aes(fill = value), colour = "white") +
    scale_fill_gradient(low = "darkgreen",high = "red",na.value = "white",
                        limits=c(exp_lengend_low, exp_lengend_high),
                        breaks=seq(exp_lengend_low,exp_lengend_high,length.out = 4),labels=scaleFUN)+
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
          axis.title.x=element_blank(),legend.text=element_text(size=7),legend.key.size = unit(0.2, "cm"),
          panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.background = element_blank(),legend.title=element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "white"),
          legend.position="right", legend.box = "vertical")
  regulators_plist[[length(regulators_plist)+1]] <- p
  multiplot(plotlist = regulators_plist,cols = 1,layout = layout)
}

#' Plot the correlation of each group and auto detected knee point. It can be used to determined which clustered are kept for further analysis.
#' 
#' @param avg_cor_list The average within group Pearson correlation coefficient of each group from gnet().
#' 
#' @return A list of indices of the data point with correlation higher than the knee point.
#' @export
plot_group_correlation <- function(gnet_result){
  avg_cor_list = gnet_result$group_score
  avg_cor_list2 <- sort(avg_cor_list,decreasing = T)

  kp <- kneepointDetection(avg_cor_list2)
  plot(1:length(avg_cor_list2),avg_cor_list2,col=c(rep(3,kp),rep(2,length(avg_cor_list2)-kp)),
       pch=1,cex =0.6,xlab='Cluster number',ylab='Average correlation',main='Cluster number vs. Average correlation')

  k1 <- avg_cor_list2[1:(kp)]
  k2 <- 1:kp
  if(kp>1){
    f1 <- lm(k1 ~ k2)
    lines(x=k2, y=predict(f1, newdata=data.frame(x=k2)),col=3,lwd=2)
  }
  l1 <- avg_cor_list2[(kp+1):length(avg_cor_list2)]
  l2 <- (kp+1):length(avg_cor_list2)
  if(length(avg_cor_list2)-kp>1){
    f2 <- lm(l1 ~ l2)
    lines(x=l2, y=predict(f2, newdata=data.frame(x=l2)),col=2,lwd=2)
  }
  return(which(avg_cor_list >= avg_cor_list2[kp]))
}

#' Build regulation modules by iteratively perform TF assigning and Gene assigning, until the assignment of genes did not change, or max number of iterations reached
#' 
#' @param input A p by n matrix of expression data of p genes and n samples, for example log2 RPKM from RNA-Seq.
#' @param reg_names A list of name of potential upstream regulators such as transcription factors.
#' @param init_group_num Initial number of function clusters used by the algorithm.
#' @param max_partition_level max_partition_level Maximum partition level in the tree.
#' @param cor_cutoff  Cutoff for within group Pearson correlation coefficient, if all data belong to a node have average correlation greater or equal to this, the node would not split anymore.
#' @param min_divide_size Minimum number of data belong to a node allowed for further split of the node.
#' @param min_divide_size Minimum number of genes allowed in a group.
#' @param max_iter Maxumum number of iterations allowed if not converged.
#' 
#' @return A list of expression data of genes, expression data of regulators, within group score, table of tree structure and final assigned group of each gene.
#' @examples
#' set.seed(1)
#' gene_num = 1000
#' tf_num = 100
#' input <- matrix(rnorm(gene_num*12),gene_num,12)
#' tf_list <- paste0('TF',1:tf_num)
#' rownames(input) <- c(tf_list,paste0('gene',1:(nrow(input)-length(tf_list))))
#' colnames(input) <- paste0('sample_',1:ncol(input))

#' gnet_result = gnet(input,tf_list,init_group_num,max_partition_level,
#' cor_cutoff,min_divide_size,min_group_size,max_iter,min_group_num)
#' @export
gnet <- function(input,reg_names,init_group_num = 20,max_partition_level = 4,
                 cor_cutoff = 0.9,min_divide_size = 3,min_group_size = 2,max_iter = 5){
  gene_data = input[!rownames(input)%in%reg_names,]
  regulator_data = input[reg_names,]
  result_all = run_gnet(gene_data,regulator_data,init_group_num,max_partition_level,cor_cutoff,min_divide_size,
                        min_group_size,max_iter)
  reg_group_table = result_all[[1]]
  gene_group_table = result_all[[2]]
  avg_cor_list <- rep(0,length(unique(gene_group_table$group)))
  for(i in 1:length(avg_cor_list)){
    cor_m <- cor(t(gene_data[gene_group_table$group==(i-1),]))
    avg_cor_list[i] <- mean(cor_m[upper.tri(cor_m)])
  }
  return(list('gene_data'=gene_data,'regulator_data'=regulator_data,'group_score' = avg_cor_list,
              'reg_group_table'=reg_group_table,'gene_group_table'=gene_group_table))
}





