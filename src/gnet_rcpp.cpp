#include <Rcpp.h>
using namespace Rcpp;

double sample_sd(NumericVector X) {
    double n = X.size();
    double sd1 = sqrt((n-1)/n)*sd(X);
    return sd1;
}

double pearsoncoeff(NumericVector X, NumericVector Y)
{
    
    return (mean(X*Y)-mean(X)*mean(Y)) / (sample_sd(X)* sample_sd(Y));
}

double rowwise_avg_cor(NumericMatrix Y,LogicalVector row_slt) {
    NumericVector cor_list((sum(row_slt)-1)*sum(row_slt)/2);
    int k =0;
    for(int i=0;i<(Y.rows()-1);++i) {
        for(int j=i+1;j<Y.rows();++j) {
            if(row_slt[i] & row_slt[j]) {
                cor_list[k] = pearsoncoeff(Y(i,_),Y(j,_));
                k++;
            }
        }
    }
    return mean(cor_list);
}

double update_score(NumericMatrix X,LogicalVector left_idx_bool,LogicalVector right_idx_bool) {
    double score_nosplit = 0;
    double score_left = 0;
    double score_right = 0;
    LogicalVector all_valid_idx_bool = left_idx_bool | right_idx_bool;
    for(int i=0;i<X.ncol();++i) {
        NumericVector x_i = X(_,i);
        NumericVector x_i_all = x_i[all_valid_idx_bool];
        score_nosplit += sum(dnorm(x_i_all,mean(x_i_all),sd(x_i_all),true));
        
        if (sum(left_idx_bool)>1){
            NumericVector x_i_left = x_i[left_idx_bool];
            score_left += sum(dnorm(x_i_left,mean(x_i_left),sd(x_i_left),true));        
        }
        if (sum(right_idx_bool)>1){
            NumericVector x_i_right= x_i[right_idx_bool];
            score_right += sum(dnorm(x_i_right,mean(x_i_right),sd(x_i_right),true)); 
        }
    }
    return score_left+score_right-score_nosplit;
}


IntegerVector find_new_split(NumericMatrix X, NumericMatrix Y,IntegerVector groups, LogicalVector feature_remaining) {
    int best_feature = -1;
    double best_score = -pow(10,6);
    IntegerVector best_split(X.nrow());
    IntegerVector groups_unique = unique(groups);
    groups_unique = groups_unique[groups_unique!=-1];
    for (IntegerVector::iterator current_branch= groups_unique.begin(); current_branch != groups_unique.end(); ++current_branch){
        LogicalVector current_split_group = groups== *current_branch;
        
        for(int feature_idx=0;feature_idx<X.ncol();++feature_idx) {
            if(!feature_remaining[feature_idx]){
                NumericVector current_feature_vals = X(_,feature_idx);
                NumericVector valid_vals = current_feature_vals[current_split_group];
                NumericVector valid_vals_unique = unique(valid_vals);
                for(NumericVector::iterator i = valid_vals_unique.begin(); i != valid_vals_unique.end(); ++i) {
                    double current_val = *i;
                    LogicalVector left_idx_bool = (current_feature_vals <= current_val) & current_split_group;
                    LogicalVector right_idx_bool = (current_feature_vals>current_val) & current_split_group;
                    double current_split_score = update_score(Y, left_idx_bool, right_idx_bool);
                    if(current_split_score > best_score) {
                        best_feature = feature_idx;
                        best_score = current_split_score;
                        best_split[left_idx_bool] = 0;
                        best_split[right_idx_bool] = 1;
                        best_split[!current_split_group] = -1;
                    }
                }
            }
        }
    }
    
    IntegerVector new_group_table_row(best_split.size()+1);
    new_group_table_row[0] = best_feature;
    new_group_table_row[Range(1, new_group_table_row.size()-1)] = best_split;
    return new_group_table_row;
}

//' Fit a regression tree.
//' 
//' Fit a regression tree based on Gaussian Likelihood score.
//' @param X A n by p matrix as input.
//' @param Y A n by q matrix as response.
//' @param max_depth Maximum depth of the tree.
//' @param cor_cutoff Cutoff for within group Pearson correlation coefficient, if all data belong to a node 
//' have average correlation greater or equal to this, the node would not split anymore.
//' @param min_divide_size Minimum number of data belong to a node allowed for further split of the node.
//' 
//' @return A matrix for sample informatrion for each partition level. First column is feature index used 
//' by the node and second is the value used to split, the rest of the columns are the split of sample: 0 means 
//' less or equal, 1 means greater and -1 means the sample does not belong to this node.
//' @examples
//' build_moduleR(X = matrix(rnorm(5*10),5,10), Y = matrix(rnorm(5*10),5,10),
//'                          max_depth=3,cor_cutoff=0.9,min_divide_size=3)
//' @export
// [[Rcpp::export]]
NumericMatrix build_module(NumericMatrix X,NumericMatrix Y, int max_depth,double cor_cutoff, int min_divide_size) {
    LogicalVector feature_remaining(X.ncol());
    int feature_num = 0;
    IntegerVector subgroup_indicator = rep(0,X.nrow());
    NumericMatrix group_table(max_depth,X.nrow()+1);
    IntegerVector groups(X.nrow());
    while (feature_num < max_depth && is_false(all(feature_remaining == TRUE)) && is_false(all(groups == -1))){
        IntegerVector new_group_table_row = find_new_split(X, Y, groups, feature_remaining);
        group_table(feature_num,_) = new_group_table_row;
        feature_num++;
        feature_remaining[new_group_table_row[0]] = TRUE;
        
        // update groups
        IntegerVector group_info = new_group_table_row[Range(1, new_group_table_row.size()-1)];
        LogicalVector left_idx = group_info==0;
        LogicalVector right_idx = group_info==1;
        int next_group_idx = max(groups)+1;
        double avg_cor_left = rowwise_avg_cor(Y,left_idx);
        double avg_cor_right = rowwise_avg_cor(Y,right_idx);
        
        if((sum(left_idx) <= min_divide_size) | (avg_cor_left > cor_cutoff)){
            groups[left_idx] = -1;
        } else {
            groups[left_idx] = next_group_idx;
            next_group_idx++;
        }
        if((sum(right_idx) <= min_divide_size) | (avg_cor_right > cor_cutoff)){
            groups[right_idx] = -1;
        } else {
            groups[right_idx] = next_group_idx;
        }
        
    }
    return group_table;
}
