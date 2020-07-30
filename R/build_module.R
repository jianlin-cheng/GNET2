#' Calculate Gaussian Likelihood score.
#' 
#' @param x A n by p matrix.
#' @param group_labels A vector of length n, indicating the group of rows.
#' @return The sum of log likelihood score of each group on each column.
#' @examples
#' calc_likelihood_score(x = matrix(rnorm(5*10),5,10), group_labels = c(rep(1,2),rep(2,3)))
#' @export
calc_likelihood_score <- function(x,group_labels){
    score_total <- 0
    for(i in unique(group_labels)){
        if(sum(group_labels==i) > 1){
            x_i <- x[group_labels==i,,drop=FALSE]
            mean_x_i <- rep(colMeans(x_i), each = nrow(x_i))
            sd_x_i = rep(colSds(x_i), each = nrow(x_i))
            density_scores <- dnorm(x_i, mean_x_i, sd_x_i, log = TRUE)
            density_scores[is.infinite(density_scores)]=0 # subgroup with all 0 in gene i are set to 0
            score_total <- score_total + sum(density_scores)
        }
    }
    return(score_total)
}

calc_likelihood_vec <- function(x,group_labels){
    score_total <- 0
    for(i in unique(group_labels)){
        if(sum(group_labels==i) > 1){
            x_i <- x[group_labels==i]
            score_total <- score_total+sum(dnorm(x_i,mean = mean(x_i),sd = sd(x_i),log = TRUE))
        }
    }
    return(score_total)
}

calc_correlation <- function(x){
    # remove any columns with 0 variance
    x <- x[,apply(x,2,var)>0]
  
    if(nrow(x)<=2 | ncol(x)<=2){
        # should not go further divide
        return(1)
    }else{
        x_cor <- cor(x)
        return(mean(abs(x_cor[upper.tri(x_cor)])))
    }
}

#' Calculate correlation within each group.
#' 
#' Calculate Pearson correlation coefficient within each group.
#' @param x A n by p matrix.
#' @param group_labels A vector of length n, indicating the group of rows.
#' @return An array of Pearson correlation coefficient for each row, rows belong to the same group have same values.
#' @examples
#' get_correlation_list(x = matrix(rnorm(5*10),5,10), group_labels = c(rep(1,2),rep(2,3)))
#' @export
get_correlation_list <- function(x,group_labels){
    cor_list <- rep(0,length(group_labels))
    for(i in unique(group_labels)){
        if(sum(group_labels==i)==1){
            cor_list[group_labels==i] <- 0
        }else{
            cor_list[group_labels==i] <- calc_correlation(x[,group_labels==i,drop=FALSE])
        }
    }
    return(cor_list)
}

get_leaf_group_labels <- function(group_table,format_plot=FALSE){
    if(!format_plot){
        group_table <- group_table[,2:ncol(group_table),drop=FALSE]
    }
    leaf_label <- rep(0,ncol(group_table))
    next_label <- 0
    for(i in seq_len(nrow(group_table))){
        current_label <- group_table[i,]
        for(j in seq_len(2)-1){
            leaf_label[current_label==j] <- next_label
            next_label <- next_label+1
        }
    }
    return(leaf_label)
}

#' Build regression tree.
#'
#' Build regression tree based on Gaussian Likelihood score.
#' @param X A n by p matrix as input.
#' @param Y A n by q matrix as response.
#' @param max_depth Maximum depth of the tree.
#' @param cor_cutoff Cutoff for within group Pearson correlation coefficient, if all data belong to 
#' a node have average correlation greater or equal to this, the node would not split anymore.
#' @param min_divide_size Minimum number of data belong to a node allowed for further split of the node.
#' 
#' @return A matrix for sample information for each tree level. First column is feature index used by the 
#' node and second is the value used to split, 
#' the rest of the columns are the split of sample: 0 means less or equal, 1 means greater and -1 
#' means the sample does not belong to this node.
#' @examples
#' build_moduleR(X = matrix(rnorm(5*10),5,10), Y = matrix(rnorm(5*10),5,10),
#'                             max_depth=3,cor_cutoff=0.9,min_divide_size=3)
#' @export
build_moduleR <- function(X,Y,max_depth,cor_cutoff,min_divide_size){
    feature_remaining <- seq_len(ncol(X))
    feature_num <- 0
    subgroup_indicator <- rep(0,nrow(X))
    groups <- c(0)
    group_table <- matrix(-1, nrow = max_depth, ncol = nrow(X)+1)
    while (feature_num < max_depth & length(groups)>0 & length(feature_remaining)>0){
        best_score <- (-(10^6))
        best_feature <- (-1)
        for(i in groups){
            current_split_group <- subgroup_indicator == i
            y_current_split <- Y[current_split_group,,drop=FALSE]
            for (j in feature_remaining){
                feature_vals <- X[current_split_group,j]
                divide_vals <- feature_vals[feature_vals != max(feature_vals)]
                score_nosplit <- calc_likelihood_score(y_current_split,rep(1,length(feature_vals)))
                for (k in divide_vals){
                    subgroup_divide <- 1-(feature_vals <= k)
                    score <- calc_likelihood_score(y_current_split,subgroup_divide) - score_nosplit
                    if (score > best_score){
                        best_score <- score
                        best_feature <- j
                        subgroup_group_labels <- rep(-1,length(current_split_group))
                        subgroup_group_labels[current_split_group] <- subgroup_divide
                    }
                }
            }
        }
        feature_num <- feature_num+1
        group_table[feature_num,] <- c(best_feature-1,subgroup_group_labels)
        subgroup_indicator_new <- get_leaf_group_labels(group_table)
        subgroup_indicator_new[subgroup_indicator==-1] <- (-1)
        subgroup_indicator <- subgroup_indicator_new
        for (i in seq_len(2)-1){
            new_divide_i <- subgroup_group_labels==i
            if (sum(new_divide_i) <= min_divide_size){
                subgroup_indicator[new_divide_i] <- (-1)
            }else if(calc_correlation(t(Y[new_divide_i,])) >= cor_cutoff){
                subgroup_indicator[new_divide_i] <- (-1)
            }
        }
        groups <- unique(subgroup_indicator[subgroup_indicator!=-1])
        feature_remaining <- feature_remaining[feature_remaining!=best_feature]
    }
    # remove unused rows
    group_table <- group_table[apply(group_table, 1, function(x)sum(x!=-1))>0,]
    return(group_table)
}

get_group_kmeans <- function(x){
  suppressWarnings({ km <- kmeans(x,centers = 3)})
  grp <- rep(0,length(x))
  o <- order(km$centers)
  grp[km$cluster==o[1]] <- -1
  grp[km$cluster==o[3]] <- 1
  return(grp)
}

#' Build split table by K-means heuristicly.
#'
#' Build split table by K-means with 3 cluster centers for each column of X
#' @param X A n by p matrix as input.
#' 
#' @return A n by p matrix with each column consists of 3 clusters: -1 for low, 0 for mid and 1 for high
#' @examples
#' split_table <- build_split_table(matrix(rnorm(5*10),5,10))
#' @export
build_split_table <- function(X){
  split_table <- matrix(0,nrow = nrow(X),ncol = ncol(X))
  for(i in seq_len(ncol(split_table))){
    split_table[,i] <- get_group_kmeans(X[,i])
  }
  return(split_table)
}


#' Build regression tree with splits are detemined by K-means heuristicly.
#'
#' Build regression tree based on Gaussian Likelihood score. The spliting of the tree is determined heuristicly by k_means.
#' @param X A n by p matrix as input.
#' @param Y A n by q matrix as response.
#' @param max_depth Maximum depth of the tree.
#' @param cor_cutoff Cutoff for within group Pearson correlation coefficient, if all data belong to 
#' a node have average correlation greater or equal to this, the node would not split anymore.
#' @param min_divide_size Minimum number of data belong to a node allowed for further split of the node.
#' @param split_table split table generated by K-means with build_split_table()
#' 
#' @return A matrix for sample informatrion for each tree level. First column is feature index used by the 
#' node and second is the value used to split, 
#' the rest of the columns are the split of sample: 0 means less or equal, 1 means greater and -1 
#' means the sample does not belong to this node.
#' @examples
#' X <- matrix(rnorm(5*10),5,10)
#' build_moduleR_heuristic(X = X, Y = matrix(rnorm(5*10),5,10),max_depth=3,cor_cutoff=0.9,
#'                         min_divide_size=3,split_table = build_split_table(X))
#' @export
build_moduleR_heuristic <- function(X,Y,max_depth,cor_cutoff,min_divide_size,split_table){
  feature_remaining <- seq_len(ncol(X))
  feature_num <- 0
  subgroup_indicator <- rep(0,nrow(X))
  groups <- c(0)
  group_table <- matrix(-1, nrow = max_depth, ncol = nrow(X)+1)
  while (feature_num < max_depth & length(groups)>0 & length(feature_remaining)>0){
    best_score <- (-(10^6))
    best_feature <- (-1)
    for(i in groups){
      current_split_group <- subgroup_indicator == i
      y_current_split <- Y[current_split_group,,drop=FALSE]
      for (j in feature_remaining){
        feature_vals <- X[current_split_group,j]
        score_nosplit <- calc_likelihood_score(y_current_split,rep(1,length(feature_vals)))
        subgroup_divide_1 <- as.numeric(split_table[current_split_group,j]>=0) #low expressed/no change
        subgroup_divide_2 <- as.numeric(split_table[current_split_group,j]>0) #high expressed/no change
        score_split1 <- score_split2 <- (-(10^6))
        if(length(unique(subgroup_divide_1))==2){
          score_split1 <- calc_likelihood_score(y_current_split,subgroup_divide_1) - score_nosplit
        }
        if(length(unique(subgroup_divide_2))==2){
          score_split2 <- calc_likelihood_score(y_current_split,subgroup_divide_2) - score_nosplit
        }
        if(score_split1 <= score_split2){
          score <- score_split2
          subgroup_divide <- subgroup_divide_2
        }else{
          score <- score_split1
          subgroup_divide <- subgroup_divide_1
        }
        if (score > best_score){
          best_score <- score
          best_feature <- j
          subgroup_group_labels <- rep(-1,length(current_split_group))
          subgroup_group_labels[current_split_group] <- subgroup_divide
        }
      }
    }
    feature_num <- feature_num+1
    group_table[feature_num,] <- c(best_feature-1,subgroup_group_labels)
    subgroup_indicator_new <- get_leaf_group_labels(group_table)
    subgroup_indicator_new[subgroup_indicator==-1] <- (-1)
    subgroup_indicator <- subgroup_indicator_new
    for (i in seq_len(2)-1){
      new_divide_i <- subgroup_group_labels==i
      if (sum(new_divide_i) <= min_divide_size){
        subgroup_indicator[new_divide_i] <- (-1)
      }else if(calc_correlation(t(Y[new_divide_i,])) >= cor_cutoff){
        subgroup_indicator[new_divide_i] <- (-1)
      }
    }
    groups <- unique(subgroup_indicator[subgroup_indicator!=-1])
    feature_remaining <- feature_remaining[feature_remaining!=best_feature]
  }
  # remove unused rows
  group_table <- group_table[apply(group_table, 1, function(x)sum(x!=-1))>0,]
  return(group_table)
}

build_moduleR_heuristic2 <- function(X,Y,max_depth,cor_cutoff,min_divide_size){
  feature_remaining <- seq_len(ncol(X))
  feature_num <- 0
  subgroup_indicator <- rep(0,nrow(X))
  groups <- c(0)
  group_table <- matrix(-1, nrow = max_depth, ncol = nrow(X)+1)
  while (feature_num < max_depth & length(groups)>0 & length(feature_remaining)>0){
    best_score <- (-(10^6))
    best_feature <- (-1)
    for(i in groups){
      current_split_group <- subgroup_indicator == i
      y_current_split <- Y[current_split_group,,drop=FALSE]
      for (j in feature_remaining){
        feature_vals <- X[current_split_group,j]
        score_nosplit <- calc_likelihood_score(y_current_split,rep(1,length(feature_vals)))
        subgroup_divide <- suppressWarnings({as.numeric(kmeans(feature_vals,centers = 2)$cluster)-1})
        if(mean(feature_vals[subgroup_divide==0]) > mean(feature_vals[subgroup_divide==1])){
          # 0 is always left and 1 is always right
          subgroup_divide <- 1- subgroup_divide
        }
        score <- calc_likelihood_score(y_current_split,subgroup_divide) - score_nosplit
        if (score > best_score){
          best_score <- score
          best_feature <- j
          subgroup_group_labels <- rep(-1,length(current_split_group))
          subgroup_group_labels[current_split_group] <- subgroup_divide
        }
      }
    }
    feature_num <- feature_num+1
    group_table[feature_num,] <- c(best_feature-1,subgroup_group_labels)
    subgroup_indicator_new <- get_leaf_group_labels(group_table)
    subgroup_indicator_new[subgroup_indicator==-1] <- (-1)
    subgroup_indicator <- subgroup_indicator_new
    for (i in seq_len(2)-1){
      new_divide_i <- subgroup_group_labels==i
      if (sum(new_divide_i) <= min_divide_size){
        subgroup_indicator[new_divide_i] <- (-1)
      }else if(calc_correlation(t(Y[new_divide_i,])) >= cor_cutoff){
        subgroup_indicator[new_divide_i] <- (-1)
      }
    }
    groups <- unique(subgroup_indicator[subgroup_indicator!=-1])
    feature_remaining <- feature_remaining[feature_remaining!=best_feature]
  }
  # remove unused rows
  group_table <- group_table[apply(group_table, 1, function(x)sum(x!=-1))>0,]
  return(group_table)
}

assign_regul <- function(regulator_data,gene_data,gene_group_table,min_group_size,
                         max_depth,cor_cutoff,min_divide_size,heuristic = TRUE,split_table = NULL){
    X <- t(regulator_data)
    group_group_labels <- unique(gene_group_table)
    group_group_labels <- group_group_labels[group_group_labels!=-1]
    reg_group_table <- matrix(-1,nrow = max(gene_group_table)+1,ncol = ncol(regulator_data))
    group_table <- matrix(-1,nrow = max_depth*(max(gene_group_table)+1),
                          ncol = ncol(regulator_data)+2)
    i <- 0
    group_table_start <- 1
    for (group_idx in group_group_labels){
        gene_idx <- gene_group_table == group_idx
        if (sum(gene_idx) >= min_group_size){
            Y <- t(gene_data[gene_idx,,drop=FALSE])
            if(heuristic){
              group_table_i <- build_moduleR_heuristic(X,Y,max_depth,cor_cutoff,min_divide_size,split_table)
            }else{
              group_table_i <- build_module(X,Y,max_depth,cor_cutoff,min_divide_size)
            }
            if(length(ncol(group_table_i))!=0){
              group_table_i <- group_table_i[apply(group_table_i, 1, function(x)(sum(x!=0)))!=0,]
              if(length(ncol(group_table_i))!=0){
                reg_group_table[i+1,] <- get_leaf_group_labels(group_table_i)
                group_table_end <- group_table_start+nrow(group_table_i)-1
                group_table[group_table_start:group_table_end,2:ncol(group_table)] <- group_table_i
                group_table[group_table_start:group_table_end,1] <- i
                i <- i+1
                group_table_start <- group_table_end+1
              }
            }
        }
    }
    # remove unused rows
    group_table <- group_table[apply(group_table, 1, function(x)sum(x!=-1))>0,]
    return(list(group_table,reg_group_table))
}

assign_gene <- function(gene_data,reg_group_table){
    gene_group_table <- rep(-1,nrow(gene_data))
    for(gene_idx in seq_len(nrow(gene_data))){
        exp_gene <- as.numeric(gene_data[gene_idx,])
        gene_group_table[gene_idx] <- which.max(
            apply(reg_group_table, 1, function(x)calc_likelihood_vec(exp_gene,x)))-1
    }
    return(gene_group_table)
}

#' Knee point detection.
#' 
#' Detect the knee point of the array.
#' @param vect A list of sorted numbers.
#' 
#' @return The index of the data point which is the knee.
#' @examples
#' kneepointDetection(sort(c(runif(10,1,3),c(runif(10,5,10))),TRUE))
#' @export
kneepointDetection <- function (vect){
    n <- length(vect)
    MinError <- 1e+08
    MinIndex <- 1
    for (i in 2:(n - 2)){
        a <- data.frame('V1'=seq_len(i), 'V2'=vect[seq_len(i)])
        l1 <- lm(a[, 2] ~ a[, 1], data = a)
        e1 <- sum(abs(l1$residuals))
        a <- data.frame('V1'=(i + 1):n, 'V2'=vect[(i + 1):n])
        l2 <- lm(a[, 2] ~ a[, 1], data = a)
        e2 <- sum(abs(l2$residuals))
        Error <- e1 + e2
        if (MinError > Error){
            MinError <- Error
            MinIndex <- i
        }
    }
    return(MinIndex)
}


assign_first_cluster <- function(gene_data,regulator_data,max_depth,init_group_num,
                                 init_method='boosting',max_group=5,nthread = 4){
    if(init_method=='boosting'){
        ipt_mat <- matrix(0,nrow = nrow(gene_data),ncol = nrow(regulator_data))
        rownames(ipt_mat) <- rownames(gene_data)
        colnames(ipt_mat) <- rownames(regulator_data)
        for(i in seq_len(nrow(gene_data))){
            dtrain <- xgb.DMatrix(data = t(regulator_data), label=gene_data[i,])
            bst <- xgb.train(data=dtrain, max.depth=max_depth, eta=0.1, nthread = nthread,nrounds =2)
            if(length(xgb.dump(model = bst, with_stats = TRUE))!=4){
              # A constant tree: no predictors were used. 
              importance_matrix <- xgb.importance(model = bst)
              ipt_mat[i,importance_matrix$Feature] <- importance_matrix$Gain
            }
        }
    }else{
        ipt_mat <- gene_data
    }
    avg_cor_list <- c()
    for(i in 2:min(init_group_num,nrow(gene_data)-1)){
        gene_group_table <- suppressWarnings({as.numeric(kmeans(ipt_mat,centers = i)$cluster)-1})
        avg_cor_list <- c(avg_cor_list,mean(get_correlation_list(t(gene_data),gene_group_table)))
    }
    if(length(avg_cor_list)>=3){
        o <- order(avg_cor_list,decreasing = TRUE)
        y <- avg_cor_list[o]
        kn <- kneepointDetection(y)
        groups_keep <- o[seq_len(max(kn,max_group))]
    }else{
        groups_keep <- seq_len(length(avg_cor_list))
    }
    gene_group_table[!(gene_group_table %in% groups_keep)] <- (-1)
    return(gene_group_table)
}

run_gnet <- function(gene_data,regulator_data,init_method = 'boosting',init_group_num = 5,max_depth = 3,
                     cor_cutoff = 0.9,min_divide_size = 3,min_group_size = 2,max_iter = 5,heuristic = TRUE,
                     max_group=5,force_split = 0.5,nthread = 4){
    message('Determining initial group number...')
    gene_group_table <- assign_first_cluster(gene_data,regulator_data,max_depth,
                                             init_group_num,init_method,max_group,nthread = nthread)
    split_table <- build_split_table(t(regulator_data))
    message('Building module networks...')
    for (i in seq_len(max_iter)) {
        message('Iteration ',i)
        assign_reg_names <- assign_regul(regulator_data,gene_data,gene_group_table,min_group_size,
                                         max_depth,cor_cutoff,min_divide_size,heuristic,split_table)
        gene_group_table_new <- assign_gene(gene_data,assign_reg_names[[2]])
        max_group_idx <- as.numeric(names(which.max(table(gene_group_table_new))))
        if(max(table(gene_group_table_new))>round(length(gene_group_table_new)*force_split)){
            gene_group_table_largest <- assign_first_cluster(gene_data[gene_group_table_new == max_group_idx,],
                                                             regulator_data,max_depth,init_group_num,init_method,max_group)
            gene_group_table_largest[gene_group_table_largest==-1] <- (-1-max(gene_group_table_new))
            gene_group_table_new[gene_group_table_new==max_group_idx] <- gene_group_table_largest+max(gene_group_table_new)
            gene_group_table_new <- as.numeric(as.factor(gene_group_table_new))
            group_to_remove <- as.numeric(names(table(gene_group_table_new))[table(gene_group_table_new) <= min_group_size])
            gene_group_table_new[gene_group_table_new %in% group_to_remove] <- -1
            gene_group_table_new <- as.numeric(as.factor(gene_group_table_new))-1
        }
        
        if(all(length(gene_group_table)==length(gene_group_table_new)) && all(gene_group_table==gene_group_table_new)){
            message('Converged.')
            break
        }else{
            gene_group_table <- gene_group_table_new
        }
    }
    message('Generating final network modules...')
    assign_reg_names <- assign_regul(regulator_data,gene_data,gene_group_table,min_group_size,
                                     max_depth,cor_cutoff,min_divide_size,heuristic,split_table)
    gene_group_table <- assign_gene(gene_data,assign_reg_names[[2]])
    message('Done.')
    tree_table_all <- assign_reg_names[[1]]
    colnames(tree_table_all) <- c('group','feature',colnames(gene_data))
    group_order <- order(gene_group_table)
    gene_list_all <- data.frame('gene'  = rownames(gene_data)[group_order],
                                'group' = gene_group_table[group_order])
    return(list('reg_group_table'=tree_table_all,'gene_group_table'=gene_list_all))
}

#' Run GNET2
#' 
#' Build regulation modules by iteratively perform regulator assigning and Gene assigning, until the 
#' assignment of genes did not change, or max number of iterations reached.
#' @param input A SummarizedExperiment object, or a p by n matrix of expression data of p genes and n samples,
#' for example log2 RPKM from RNA-Seq.
#' @param reg_names A list of potential upstream regulators names, for example a list of known transcription factors.
#' @param init_method Cluster initialization, can be "boosting" or "kmeans", default is using "boosting".
#' @param init_group_num Initial number of function clusters used by the algorithm.
#' @param max_depth max_depth Maximum depth of the tree.
#' @param cor_cutoff    Cutoff for within group Pearson correlation coefficient, if all data belong to a node have
#' average correlation greater or equal to this, the node would not split anymore.
#' @param min_divide_size Minimum number of data belong to a node allowed for further split of the node.
#' @param min_group_size Minimum number of genes allowed in a group.
#' @param max_iter Maxumum number of iterations allowed if not converged.
#' @param heuristic If the splites of the regression tree is determined by k-means heuristicly.
#' @param max_group Max number of group allowed for the first clustering step, default equals init_group_num and is set to 0.
#' @param force_split Force split the largest gene group into smaller groups by kmeans. Default is 0.5(Split if it contains more than half target genes)
#' @param nthread Number of threads to run GBDT based clustering
#' 
#' @return A list of expression data of genes, expression data of regulators, within group score, table of tree 
#' structure and final assigned group of each gene.
#' @examples
#' set.seed(1)
#' init_group_num = 8
#' init_method = 'boosting'
#' exp_data <- matrix(rnorm(50*10),50,10)
#' reg_names <- paste0('TF',1:5)
#' rownames(exp_data) <- c(reg_names,paste0('gene',1:(nrow(exp_data)-length(reg_names))))
#' colnames(exp_data) <- paste0('condition_',1:ncol(exp_data))
#' se <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=exp_data))
#' gnet_result <- gnet(se,reg_names,init_method,init_group_num)
#' @export
gnet <- function(input,reg_names,init_method= 'boosting',init_group_num = 4,max_depth = 3,
                 cor_cutoff = 0.9,min_divide_size = 3,min_group_size = 2,max_iter = 5,
                 heuristic = TRUE,max_group = 0,force_split = 0.5,nthread = 4){
    if(is(input,class2 = "SummarizedExperiment")){
        input <- assay(input)
    }
    if(max_group == 0){
      max_group <- init_group_num
    }
    input <- input[apply(input, 1, var) > 0.0001,]
    gene_data <- input[!rownames(input)%in%reg_names,,drop=FALSE]
    regulator_data <- input[reg_names,,drop=FALSE]
    result_all <- run_gnet(gene_data,regulator_data,init_method,init_group_num,max_depth,cor_cutoff,
                           min_divide_size,min_group_size,max_iter,heuristic,max_group,force_split,nthread = nthread)
    reg_group_table <- result_all[[1]]
    gene_group_table <- result_all[[2]]
    
    # sanity check: remove all modules without any genes assigned
    
    
    
    group_idx <- unique(reg_group_table[,1])
    reg_group_table_filtered <- gene_group_table_filtered <- NULL
    avg_cor_list <- c()
    current_group_idx <- 1
    regulators <- target_genes <- list()
    for (i in seq_len(length(group_idx))) {
        if(sum(gene_group_table[,2]==group_idx[i])>=min_group_size){
            current_tree <- reg_group_table[reg_group_table[,1] == group_idx[i],]
            current_tree[,1] <- current_group_idx
            current_gene_group <- gene_group_table[gene_group_table$group==group_idx[i],]
            current_gene_group$group<- current_group_idx
            reg_group_table_filtered <- rbind(reg_group_table_filtered,current_tree)
            gene_group_table_filtered <- rbind(gene_group_table_filtered,current_gene_group)
            
            cor_m <- cor(t(gene_data[current_gene_group$gene,,drop=FALSE]))
            avg_cor_list <- c(avg_cor_list,mean(cor_m[upper.tri(cor_m)]))
            
            regulators[[i]] <- rownames(regulator_data)[current_tree[,2]+1]
            target_genes[[i]] <- current_gene_group$gene
            current_group_idx <- current_group_idx + 1
        }
    }
    
    if(nrow(reg_group_table_filtered)<=2)warning('Too few modules generated, you may wish to try with higher cor_cutoff.')
    return(list('gene_data' = gene_data,'regulator_data' = regulator_data,'group_score' = avg_cor_list,
                'reg_group_table' = reg_group_table_filtered,'gene_group_table' = gene_group_table_filtered,
                'modules_count' = current_group_idx-1,'regulators' = regulators,
                'target_genes' = target_genes))
}

sum_scores <- function(el_all,el_input){
  scores <- rep(0,nrow(el_all))
  for (i in seq_len(nrow(el_all))) {
    idx <- (el_input[,1] == el_all[i,1] & el_input[,2] == el_all[i,2]) | (el_input[,2] == el_all[i,1] & el_input[,1] == el_all[i,2])
    scores[i] <- sum(abs(el_input[idx,3]))
  }
  return(scores)
}

#' Extract the network from the gnet result
#' 
#' Extract the network as edge list from the gnet result. For a module, each regulator and downstream gene will form a directed edge.
#' @param gnet_result Returned results from gnet().
#' 
#' @return A matrix of scores of for the regulator-target interaction.
#' @examples
#' set.seed(1)
#' init_group_num = 8
#' init_method = 'kmeans'
#' exp_data <- matrix(rnorm(50*10),50,10)
#' reg_names <- paste0('TF',1:5)
#' rownames(exp_data) <- c(reg_names,paste0('gene',1:(nrow(exp_data)-length(reg_names))))
#' colnames(exp_data) <- paste0('condition_',1:ncol(exp_data))
#' se <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=exp_data))
#' gnet_result <- gnet(se,reg_names,init_method,init_group_num)
#' edge_list <- extract_edges(gnet_result)
#' @export
extract_edges <- function(gnet_result){
  n_el <- sum(sapply(gnet_result$regulators, length) *sapply(gnet_result$target_genes, length))+ 
    sum(sapply(gnet_result$regulators, length) *sapply(gnet_result$regulators, length))
  el <- data.frame(matrix(0,nrow = n_el,ncol = 3))
  current_start = 1
  for (i in seq_len(gnet_result$modules_count)) {
    nr <- length(gnet_result$regulators[[i]])
    el[current_start:(current_start+nr^2-1),1:2] <- expand.grid(gnet_result$regulators[[i]],gnet_result$regulators[[i]],stringsAsFactors =FALSE)
    el[current_start:(current_start+nr^2-1),3] <- gnet_result$group_score[i]
    current_start <- current_start+nr^2
    
    nr <- length(gnet_result$regulators[[i]]) * length(gnet_result$target_genes[[i]])
    el[current_start:(current_start+nr-1),1:2] <- expand.grid(gnet_result$regulators[[i]],gnet_result$target_genes[[i]],stringsAsFactors =FALSE)
    el[current_start:(current_start+nr-1),3] <- gnet_result$group_score[i]
    current_start <- current_start+nr
  }
  
  g_all <- graph_from_edgelist(as.matrix(el[,c(1,2)]),directed = FALSE)
  el1a <- unique(get.edgelist(g_all))
  el1b <- sum_scores(el1a,el)
  el <- cbind.data.frame(el1a,'score'=el1b)
  
  colnames(el)<- c('regulator','target','score')
  el1 <- el %>% group_by_(.dots = c('regulator','target')) %>% summarise_all(list('score' = sum))
  el2 <- data.frame('regulator'=as.character(el1$regulator),
                    'target'=as.character(el1$target),
                    'score'=as.numeric(el1$score),stringsAsFactors = FALSE)
  el2 <- el2[el1$regulator!= el2$target,]
  el2$score <- el2$score/max(el2$score)
  rownames(el2) <- seq_len(nrow(el2))
  
  
  result_gnet <- dcast(el,regulator~target,value.var="score",fill = 0)
  rownames(result_gnet) <- result_gnet$regulator
  result_gnet <- result_gnet[,2:ncol(result_gnet)]
  
  reg_list <- rownames(gnet_result$regulator_data)
  target_list <- c(reg_list,rownames(gnet_result$gene_data))
  
  mat  <- matrix(0,nrow = length(reg_list),ncol = length(target_list))
  rownames(mat) <- reg_list
  colnames(mat) <- target_list
  
  
  cor_all <- cor(t(rbind(gnet_result$regulator_data,gnet_result$gene_data)),
                 method = 'spearman', use = "complete.obs")
  diag(cor_all) <- 0
  cor_all2 <- cor_all[rownames(result_gnet),colnames(result_gnet)]
  cor_all2[is.na(cor_all2)] <- 0
  cor_all2 <- cor_all2^2
  cor_all2[cor_all2>1-1e-6] <- 1-1e-6
  
  v <- as.matrix(result_gnet - 2 * log(1 - cor_all2))
  
  
  
  mat[rownames(v),colnames(v)] <- v
  return(mat)
  
}

get_sample_cluster <- function(group_table,group_idx){
  group_table <- group_table[group_table[,1]==group_idx,3:ncol(group_table),drop=FALSE]
  leaf_label <- rep(0,ncol(group_table))
  next_label <- 0
  for(i in seq_len(nrow(group_table))){
    current_label <- group_table[i,]
    for(j in c(0,1)){
      leaf_label[current_label==j] <- next_label
      next_label <- next_label+1
    }
  }
  return(as.numeric(factor(leaf_label)))
}

ari <- function (x, y){
  tab <- table(x, y)
  if (all(dim(tab) == c(1, 1))) 
    return(1)
  a <- sum(choose(tab, 2))
  b <- sum(choose(rowSums(tab), 2)) - a
  c <- sum(choose(colSums(tab), 2)) - a
  d <- choose(sum(tab), 2) - a - b - c
  ARI <- (a - (a + b) * (a + c)/(a + b + c + d))/((a + b + 
                                                     a + c)/2 - (a + b) * (a + c)/(a + b + c + d))
  return(ARI)
}

similarity_score_unranked <- function(gnet_result,group){
  s_list<- rep(0,gnet_result$modules_count)
  for (i in seq_len(gnet_result$modules_count)) {
    s_list[i] <- ari(get_sample_cluster(gnet_result$reg_group_table,i),group)
  }
  return(s_list)
}

similarity_score_ranked <- function(gnet_result,group){
  kld <- rep(0,gnet_result$modules_count)
  for (i in seq_len(gnet_result$modules_count)) {
    group_module <- get_sample_cluster(gnet_result$reg_group_table,i)
    group_module_dist <- as.matrix(dist(group_module))
    group_dist <- as.matrix(dist(group))
    
    p <- dnorm(group_dist*group_dist,0,sqrt(max(group_dist)))
    q <- dnorm(group_module_dist*group_module_dist,0,sqrt(max(group_module_dist)))
    
    
    for (j in seq_len(nrow(group_dist))) {
      p[j,] <- p[j,]/(sum(p[j,])-p[j,])
      q[j,] <- q[j,]/(sum(q[j,])-q[j,])
    }
    
    p1 <- p[upper.tri(p)]
    q1 <- q[upper.tri(q)]
    
    kld[i] <- sum(p1*log(p1/q1))
  }
  
  s_list <- -(kld-mean(kld))/sd(kld)
  return(s_list)
}

#' Compute the similarity from a predefined condition group
#' 
#' Compute the similarity between a predefined condition grouping and the sample cluster of each module, which is defined
#' as the Adjusted Rand index between the two vectors, or the inverse of K-L divergence between the upper triangle matrix 
#' of the pairwise distance of predefined ranked condition grouping and the pairwise distance of sample cluster of each module.
#' @param gnet_result Returned results from gnet().
#' @param group predefined condition grouping
#' @param ranked the grouping information is categorical(treatment/control) or ordinal(dosage, time points)?
#' 
#' @return A list of similarity scores between a predefined condition grouping and the sample cluster of each module.
#' @examples
#' set.seed(1)
#' init_group_num = 8
#' init_method = 'kmeans'
#' exp_data <- matrix(rnorm(50*10),50,10)
#' reg_names <- paste0('TF',1:5)
#' rownames(exp_data) <- c(reg_names,paste0('gene',1:(nrow(exp_data)-length(reg_names))))
#' colnames(exp_data) <- paste0('condition_',1:ncol(exp_data))
#' se <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=exp_data))
#' gnet_result <- gnet(se,reg_names,init_method,init_group_num)
#' s <- similarity_score(gnet_result,rep(1:5,each = 2))
#' @export
similarity_score <- function(gnet_result,group,ranked=FALSE){
    if(ranked){
      return(similarity_score_ranked(gnet_result,group))
    }else{
      return(similarity_score_unranked(gnet_result,group))
    }
}



