multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    # From http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
    plots <- c(list(...), plotlist)
    numPlots <- length(plots)
    if (is.null(layout)) {
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                                         ncol = cols, nrow = ceiling(numPlots/cols))
    }
    if (numPlots==1) {
        print(plots[[1]])
    } else {
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

        # Make each plot, in the correct location
        for (i in seq_len(numPlots)) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}

#' Plot the regression tree.
#' 
#' Plot the regression tree given the index of a module.
#' @param gnet_result Results returned by gnet().
#' @param group_idx Index of the module.
#' 
#' @return None
#' @examples
#' set.seed(1)
#' init_group_num = 5
#' init_method = 'boosting'
#' exp_data <- matrix(rnorm(50*10),50,10)
#' reg_names <- paste0('TF',1:5)
#' rownames(exp_data) <- c(reg_names,paste0('gene',1:(nrow(exp_data)-length(reg_names))))
#' colnames(exp_data) <- paste0('condition_',1:ncol(exp_data))
#' se <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=exp_data))
#' gnet_result <- gnet(se,reg_names,init_method,init_group_num)
#' plot_tree(gnet_result,group_idx=1)
#' @export
plot_tree <- function(gnet_result,group_idx){
    regulator_data <- gnet_result$regulator_data
    gene_data <- gnet_result$gene_data
    reg_group_table <- gnet_result$reg_group_table
    gene_group_table <- gnet_result$gene_group_table
    tree_i <- reg_group_table[reg_group_table[,1]==group_idx,,drop=FALSE]
    label_list <- c()
    for(i in seq_len(nrow(tree_i))){
        coverage_i <- sum(tree_i[i,3:ncol(tree_i)]>=0)
        feature_name_i <- rownames(regulator_data)[tree_i[i,2]+1]
        split_i <- round(max(regulator_data[feature_name_i,tree_i[i,3:ncol(tree_i)]==0]),4)
        label_list <- c(label_list,paste0(feature_name_i,'\n <= ',split_i,'\nCoverage: ',coverage_i))
    }
    leaf_idx <- 0
    from_list <- to_list <- edge_labels <- c()
    for(i in seq_len((nrow(tree_i)-1))){
        for(j in seq_len(2)-1){
            found_node <- FALSE
            split_i_j <- tree_i[i,3:ncol(tree_i)]==j
            for (k in (i+1):nrow(tree_i)) {
                node_k <- tree_i[k,3:ncol(tree_i)]>=0
                if(identical(split_i_j,node_k)){
                    from_list <- c(from_list,i)
                    to_list <- c(to_list,k)
                    edge_labels <- c(edge_labels,as.character(j==0))
                    found_node <- TRUE
                }
            }
            if(!found_node){
                label_list <- c(label_list,paste('leaf',leaf_idx))
                from_list <- c(from_list,i)
                to_list <- c(to_list,length(label_list))
                edge_labels <- c(edge_labels,as.character(j==0))
                leaf_idx <- leaf_idx+1
            }
        }
    }
    label_list <- c(label_list,paste('leaf',leaf_idx),paste('leaf',leaf_idx+1))
    from_list <- c(from_list,nrow(tree_i),nrow(tree_i))
    to_list <- c(to_list,length(label_list)-1,length(label_list))
    edge_labels <- c(edge_labels,'true','false')
    node_df <- create_node_df(n = length(label_list),type = "a",label = label_list,style = "filled",
                                                        color = "aqua",shape = "ellipse",width = 0.9,fontsize = 16)
    edge_df <- create_edge_df(from = from_list,to = to_list,label = tolower(edge_labels),fontsize = 16)
    graph <- create_graph(nodes_df = node_df,edges_df = edge_df,attr_theme = NULL) 
    # graph <- add_global_graph_attrs(graph,attr_type = "graph",attr = c("layout", "rankdir","fontsize"),
    #                                                                  value = c("dot", "LR"))
    # graph <- add_global_graph_attrs(graph,attr_type = "node", attr = c("fillcolor", "style", "fontname"), 
    #                                         value = c("Azure","filled", "Helvetica"))
    render_graph(graph)
}

get_row_order <- function(group_table1,regulator_data1,exp_data1){
    for(i in seq_len(nrow(group_table1))){
      row_order <- seq_len(ncol(group_table1))
      group_info <- group_table1[i,]
      current_parent <- group_info %in% c(0,1)
      row_order[current_parent] <- c(row_order[group_info==0],row_order[group_info==1])
      group_table1 <- group_table1[,row_order]
      regulator_data1 <- regulator_data1[,row_order]
      exp_data1 <- exp_data1[,row_order]
    }
    return(list(group_table1,regulator_data1,exp_data1))
}

#' Plot a module
#' 
#' Plot the regulators module and heatmap of the expression inferred downstream genes for each sample. 
#' It can be interpreted as two parts: the bars at the top shows how samples are splited by the 
#' regression tree and the heatmap at the bottom shows how downstream genes are regulated by each 
#' subgroup determined by the regulators.
#' @param gnet_result Results returned by gnet().
#' @param group_idx Index of the module.
#' @return None
#' @examples
#' set.seed(1)
#' init_group_num = 5
#' init_method = 'boosting'
#' exp_data <- matrix(rnorm(50*10),50,10)
#' reg_names <- paste0('TF',1:5)
#' rownames(exp_data) <- c(reg_names,paste0('gene',1:(nrow(exp_data)-length(reg_names))))
#' colnames(exp_data) <- paste0('condition_',1:ncol(exp_data))
#' se <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=exp_data))
#' gnet_result <- gnet(se,reg_names,init_method,init_group_num)
#' plot_gene_group(gnet_result,group_idx=1)
#' @export
plot_gene_group <- function(gnet_result,group_idx){
    gene_data <- gnet_result$gene_data
    regulator_data <- gnet_result$regulator_data
    reg_group_table <- gnet_result$reg_group_table
    gene_group_table <- gnet_result$gene_group_table
    exp_data1 <- gene_data[gene_group_table$gene[gene_group_table$group==group_idx],,drop=FALSE]
    regulator_data1 <- regulator_data[reg_group_table[reg_group_table[,1]==group_idx,2]+1,]
    group_table1 <- reg_group_table[reg_group_table[,1]==group_idx,3:ncol(reg_group_table)]
    leaf_labels <- get_leaf_group_labels(group_table1,format_plot = TRUE)
    
    row_order_list <- get_row_order(group_table1,regulator_data1,exp_data1)
    group_table2 <- row_order_list[[1]]
    regulator_data2 <- row_order_list[[2]]
    exp_data2 <- row_order_list[[3]]
    
    # group_table2 <- group_table1[,row_order,drop=FALSE]
    # regulator_data2 <- regulator_data1[,row_order,drop=FALSE]
    # exp_data2 <- exp_data1[,row_order,drop=FALSE]
    test_regulators_names <- rownames(regulator_data2)
    layout=matrix(c(seq_len(length(test_regulators_names)),rep(length(test_regulators_names)+1,
                                                               length(test_regulators_names)*2)),ncol=1)
    regulators_plist <- list()
    scaleFUN <- function(x) sprintf("%.3f", x)

    # add Regulators bars
    for(i in seq_len(length(test_regulators_names))){
        reg_data_mask <- group_table2[i,]==-1
        exp_val <- as.numeric(regulator_data2[i,])
        exp_val[reg_data_mask] <- NA
        lengend_low <- min(exp_val,na.rm = TRUE)
        lengend_high <- max(exp_val,na.rm = TRUE)
        exp_val1 <- rbind.data.frame(matrix(NA,nrow = 1,ncol = length(exp_val)),exp_val,stringsAsFactors=FALSE)
        rownames(exp_val1) <- seq_len(nrow(exp_val1))
        exp_val.m <- melt(exp_val1,id.vars = NULL)
        exp_val.m <- cbind.data.frame('y_idx'=rep(seq_len(nrow(exp_val1)),ncol(exp_val1)),
                                      exp_val.m,stringsAsFactors=FALSE)
        exp_label <- rep('',ncol(exp_val1))
        if(which(group_table2[i,]==0)[1] < which(group_table2[i,]==1)[1]){
          # left low, right high
          exp_label[max(which(group_table2[i,]==0))] <- '<- | Low'
          exp_label[min(which(group_table2[i,]==1))] <- 'High | ->'
        }else{
          # left high, right low
          exp_label[max(which(group_table2[i,]==1))] <- '<- | High'
          exp_label[min(which(group_table2[i,]==0))] <- 'Low | ->'
        }
        p <- ggplot(exp_val.m, aes_string('variable', 'y_idx')) + 
            geom_tile(aes_string(fill = 'value'), colour = "white") +
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
    exp_data2 <- as.data.frame(t(scale(t(exp_data2))))
    suppressWarnings({d <- dist(exp_data2, method = "euclidean")})
    fit <- hclust(d, method="ward")
    exp_data2 <- exp_data2[fit$order,]
    exp_lengend_low <- min(exp_data2)
    exp_lengend_high <- max(exp_data2)
    test_data.m <- melt(cbind.data.frame('gene'=rownames(exp_data2),exp_data2,stringsAsFactors=FALSE),id.vars = 'gene')
    p <- ggplot(test_data.m, aes_string('variable', 'gene')) + 
        geom_tile(aes_string(fill = 'value'), colour = "white") +
        scale_fill_gradient(low = "darkgreen",high = "red",na.value = "white",
                                                limits=c(exp_lengend_low, exp_lengend_high),
                                                breaks=seq(exp_lengend_low,exp_lengend_high,length.out = 4),labels=scaleFUN)+
        theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
                    axis.title.x=element_blank(),legend.text=element_text(size=7),
                    legend.key.size = unit(0.2, "cm"),
                    panel.border = element_blank(), panel.grid.major = element_blank(),
                    panel.background = element_blank(),legend.title=element_blank(),
                    axis.text.x = element_text(angle = 45, hjust = 1),
                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "white"),
                    legend.position="right", legend.box = "vertical")
    regulators_plist[[length(regulators_plist)+1]] <- p
    multiplot(plotlist = regulators_plist,cols = 1,layout = layout)
}

#' Plot the correlation of each group
#' 
#' Plot the correlation of each group and auto detected knee point. It can be used to 
#' determined which clustered are kept for further analysis.
#' @param gnet_result Results returned by gnet().
#' 
#' @return A list of indices of the data point with correlation higher than the knee point.
#' @examples
#' set.seed(1)
#' gnet_result <- list('group_score'=c(runif(10,1,3),c(runif(10,5,3))))
#' group_keep <- plot_group_correlation(gnet_result)
#' @export
plot_group_correlation <- function(gnet_result){
    avg_cor_list <- gnet_result$group_score
    avg_cor_list2 <- sort(avg_cor_list,decreasing = TRUE)
    kp <- kneepointDetection(avg_cor_list2)
    plot(seq_len(length(avg_cor_list2)),avg_cor_list2,
         col=c(rep(3,kp),rep(2,length(avg_cor_list2)-kp)),
         pch=1,cex =0.6,xlab='Cluster number',ylab='Average correlation',
         main='Cluster number vs. Average correlation')
    k1 <- avg_cor_list2[seq_len(kp)]
    k2 <- seq_len(kp)
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