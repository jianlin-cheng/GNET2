library(Rcpp)
library(microbenchmark)

sourceCpp('~/Dropbox/MU/workspace/new/src/test/c_test.cpp')

cppFunction('
int test(LogicalVector x){
    return which_max(x);
}

')


update_scoreR = function(Y,left_idx_bool,right_idx_bool){
  score_nosplit = 0
  score_left = 0
  score_right = 0
  all_valid_idx_bool = left_idx_bool | right_idx_bool
  for(i in 1:ncol(Y)){
    y_i = Y[,i]
    y_i_all = y_i[all_valid_idx_bool]
    y_i_left = y_i[left_idx_bool]
    y_i_right = y_i[right_idx_bool]
    score_nosplit = score_nosplit+sum(dnorm(y_i_all,mean(y_i_all),sd(y_i_all),T))
    score_left = score_left+sum(dnorm(y_i_left,mean(y_i_left),sd(y_i_left),T))
    score_right = score_right+sum(dnorm(y_i_right,mean(y_i_right),sd(y_i_right),T))
  }
  return(score_left+score_right-score_nosplit)
}

find_new_splitR = function(X,Y,groups,feature_remaining){
  best_score = -(10^6)
  best_split = rep(0,nrow(X))
  best_feature = -1
  groups_unique = unique(groups)
  groups_unique = groups_unique[groups_unique!=-1]
  for (current_branch in groups_unique){
    current_split_group = groups==current_branch
    for (i in 1:ncol(X)) {
      if(!feature_remaining[i]){
        # print(paste(i,best_score))
        valid_val = X[current_split_group,i]
        
        for(j in valid_val){
          left_idx_bool = (valid_val<=j)&current_split_group
          right_idx_bool = (valid_val>j)&current_split_group
          current_score = update_score(Y,left_idx_bool,right_idx_bool)
          if(current_score>best_score){
            best_score = current_score
            best_split[!current_split_group] = -1
            best_split[left_idx_bool] = 0
            best_split[right_idx_bool] = 1
            best_feature = i
          }
        }
        # print(paste(a,collapse = '; '))
      }
    }
  }
  
  return(c(best_feature-1,best_split))
}

rowwise_avg_corR = function(Y){
  cor_list =c()
  for(i in 1:(nrow(Y)-1)){
    for(j in (i+1):nrow(Y)){
      cor_list <- c(cor_list,cor(Y[i,],Y[j,]))
    }
    
  }
  return(mean(cor_list))
}


max_partition_level = 4
cor_cutoff = 0.9
min_divide_size = 3
min_group_size = 3
max_iter = 5
min_group_num = 3
init_group_num = 5

set.seed(1)

rnaseq_data <- read.csv('~/Dropbox/MU/workspace/new/data/test/rnaseq_data.csv',row.names = 1)
tf_list <- read.csv('~/Dropbox/MU/workspace/new/data/test/tf_list.csv',header = F)$V1

gene_data = rnaseq_data[!rownames(rnaseq_data)%in%tf_list,]
tf_data = rnaseq_data[tf_list,]


tf_list = rownames(tf_data)
gene_group_table = as.numeric(kmeans(gene_data,centers = init_group_num)$cluster)-1
gene_idx = gene_group_table == 2
X = t(tf_data)
group_labels = 0:max(gene_group_table)
tf_group_table = matrix(0,nrow = 0,ncol = ncol(tf_data))
Y = t(gene_data[gene_idx,])


groups = rep(0,nrow(X))
feature_remaining = rep(F,ncol(X))
aaa=find_new_split(X, Y, groups,feature_remaining)

bbb=find_new_splitR(X, Y, groups,feature_remaining)

sourceCpp('~/Dropbox/MU/workspace/new/src/test/c_test.cpp')


aaa=test1(X, Y, groups,feature_remaining)
bbb=test2(X, Y, groups,feature_remaining)


test2 = function(X,Y,groups,feature_remaining){
  best_score = -(10^6)
  best_split = rep(0,nrow(X))
  best_feature = -1
  groups_unique = unique(groups)
  groups_unique = groups_unique[groups_unique!=-1]
  current_branch = 0

  current_split_group = groups==current_branch
  i =1 
  valid_vals = X[current_split_group,i]
  j = valid_vals[1]
  left_idx_bool = (valid_val<=j)&current_split_group
  right_idx_bool = (valid_val>j)&current_split_group
  current_score = update_score(Y,left_idx_bool,right_idx_bool)
  print(current_score)
  
  return(left_idx_bool)
}
