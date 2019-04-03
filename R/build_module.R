#' Calculate Gaussian Likelihood score.
#' 
#' @param x A n by p matrix.
#' @param labels A vector of length n, indicating the group of rows.
#' @return The sum of log likelihood score of each group on each column.
#' @examples
#' calc_likelihood_score(x = matrix(rnorm(50*100),50,100), labels = c(rep(1,20),rep(2,30)))
#' @export
calc_likelihood_score <- function(x,labels){
    score_total <- 0
    for(i in unique(labels)){
        if(sum(labels==i) > 1){
            x_i <- x[labels==i,,drop=FALSE]
            mean_x_i <- rep(colMeans(x_i), each = nrow(x_i))
            sd_x_i = rep(colSds(x_i), each = nrow(x_i))
            density_scores <- dnorm(x_i, mean_x_i, sd_x_i, log = TRUE)
            score_total <- score_total + sum(density_scores)
        }
    }
    return(score_total)
}

calc_likelihood_vec <- function(x,labels){
    score_total <- 0
    for(i in unique(labels)){
        if(sum(labels==i) > 1){
            x_i <- x[labels==i]
            score_total <- score_total+sum(dnorm(x_i,mean = mean(x_i),sd = sd(x_i),log = TRUE))
        }
    }
    return(score_total)
}

calc_correlation <- function(x){
    if(nrow(x)<2 | ncol(x)<2){
        return(0)
    }else{
        x_cor <- cor(x)
        return(mean(abs(x_cor[upper.tri(x_cor)])))
    }
}

#' Calculate correlation within each group.
#' 
#' Calculate Pearson correlation coefficient within each group.
#' @param x A n by p matrix.
#' @param labels A vector of length n, indicating the group of rows.
#' @return An array of Pearson correlation coefficient for each row, rows belong to the same group have same values.
#' @examples
#' get_correlation_list(x = matrix(rnorm(50*100),50,100), labels = c(rep(1,20),rep(2,30)))
#' @export
get_correlation_list <- function(x,labels){
    cor_list <- rep(0,length(labels))
    for(i in unique(labels)){
        if(sum(labels==i)==1){
            cor_list[labels==i] <- 0
        }else{
            cor_list[labels==i] <- calc_correlation(x[,labels==i,drop=FALSE])
        }
    }
    return(cor_list)
}

get_leaf_labels <- function(group_table,format_plot=FALSE){
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
#' @return A matrix for sample informatrion for each tree level. First column is feature index used by the 
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
                        subgroup_labels <- rep(-1,length(current_split_group))
                        subgroup_labels[current_split_group] <- subgroup_divide
                    }
                }
            }
        }
        feature_num <- feature_num+1
        group_table[feature_num,] <- c(best_feature-1,subgroup_labels)
        subgroup_indicator_new <- get_leaf_labels(group_table)
        subgroup_indicator_new[subgroup_indicator==-1] <- (-1)
        subgroup_indicator <- subgroup_indicator_new
        for (i in seq_len(2)-1){
            new_divide_i <- subgroup_labels==i
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
                         max_depth,cor_cutoff,min_divide_size){
    X <- t(regulator_data)
    group_labels <- unique(gene_group_table)
    group_labels <- group_labels[group_labels!=-1]
    reg_group_table <- matrix(-1,nrow = max(gene_group_table)+1,ncol = ncol(regulator_data))
    group_table <- matrix(-1,nrow = max_depth*(max(gene_group_table)+1),
                          ncol = ncol(regulator_data)+2)
    i <- 0
    group_table_start <- 1
    for (group_idx in group_labels){
        gene_idx <- gene_group_table == group_idx
        if (sum(gene_idx) >= min_group_size){
            Y <- t(gene_data[gene_idx,,drop=FALSE])
            group_table_i <- build_moduleR(X,Y,max_depth,cor_cutoff,min_divide_size)
            group_table_i <- group_table_i[apply(group_table_i, 1, function(x)(sum(x!=0)))!=0,]
            reg_group_table[i+1,] <- get_leaf_labels(group_table_i)
            group_table_end <- group_table_start+nrow(group_table_i)-1
            group_table[group_table_start:group_table_end,2:ncol(group_table)] <- group_table_i
            group_table[group_table_start:group_table_end,1] <- i
            i <- i+1
            group_table_start <- group_table_end+1
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
#' kneepointDetection(sort(c(runif(20,1,3),c(runif(20,5,10))),TRUE))
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


assign_first_cluster <- function(gene_data,regulator_data,max_depth,
                                 init_group_num,init_method='boosting'){
    if(init_method=='boosting'){
        ipt_mat <- matrix(0,nrow = nrow(gene_data),ncol = nrow(regulator_data))
        rownames(ipt_mat) <- rownames(gene_data)
        colnames(ipt_mat) <- rownames(regulator_data)
        for(i in seq_len(nrow(gene_data))){
            dtrain <- xgb.DMatrix(data = t(regulator_data), label=gene_data[i,])
            bst <- xgb.train(data=dtrain, max.depth=max_depth, eta=0.1, nthread = 2,nrounds =2)
            importance_matrix <- xgb.importance(model = bst)
            ipt_mat[i,importance_matrix$Feature] <- importance_matrix$Gain
        }
    }else{
        ipt_mat <- gene_data
    }
    avg_cor_list <- c()
    for(i in 2:min(init_group_num,nrow(gene_data)-1)){
        gene_group_table <- as.numeric(kmeans(ipt_mat,centers = i)$cluster)-1
        avg_cor_list <- c(avg_cor_list,mean(get_correlation_list(t(gene_data),gene_group_table)))
    }
    if(length(avg_cor_list)>=3){
        o <- order(avg_cor_list,decreasing = TRUE)
        y <- avg_cor_list[o]
        kn <- kneepointDetection(y)
        groups_keep <- o[seq_len(max(kn,3))]
    }else{
        groups_keep <- seq_len(length(avg_cor_list))
    }
    gene_group_table[!(gene_group_table %in% groups_keep)] <- (-1)
    return(gene_group_table)
}

run_gnet <- function(gene_data,regulator_data,init_method = 'boosting',init_group_num = 5,
                     max_depth = 3,cor_cutoff = 0.9,min_divide_size = 3,min_group_size = 2,max_iter = 5){
    message('Determining initial group number...')
    gene_group_table <- assign_first_cluster(gene_data,regulator_data,max_depth,
                                             init_group_num,init_method)
    message('Building module networks...')
    for (i in seq_len(max_iter)) {
        message('Iteration ',i)
        assign_reg_names <- assign_regul(regulator_data,gene_data,gene_group_table,
                                         min_group_size,max_depth,cor_cutoff,min_divide_size)
        gene_group_table_new <- assign_gene(gene_data,assign_reg_names[[2]])
        if(all(length(gene_group_table)==length(gene_group_table_new)) && all(gene_group_table==gene_group_table_new)){
            break
        }else{
            gene_group_table <- gene_group_table_new
        }
    }
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
                 cor_cutoff = 0.9,min_divide_size = 3,min_group_size = 2,max_iter = 5){
    if(is(input,class2 = "SummarizedExperiment")){
        input <- assay(input)
    }
    gene_data <- input[!rownames(input)%in%reg_names,,drop=FALSE]
    regulator_data <- input[reg_names,,drop=FALSE]
    result_all <- run_gnet(gene_data,regulator_data,init_method,init_group_num,max_depth,
                                                 cor_cutoff,min_divide_size,min_group_size,max_iter)
    reg_group_table <- result_all[[1]]
    gene_group_table <- result_all[[2]]
    avg_cor_list <- rep(0,length(unique(gene_group_table$group)))
    for(i in seq_len(length(avg_cor_list))){
        cor_m <- cor(t(gene_data[gene_group_table$group==(i-1),,drop=FALSE]))
        avg_cor_list[i] <- mean(cor_m[upper.tri(cor_m)])
    }
    return(list('gene_data' = gene_data,'regulator_data' = regulator_data,'group_score' = avg_cor_list,
                            'reg_group_table' = reg_group_table,'gene_group_table' = gene_group_table))
}
