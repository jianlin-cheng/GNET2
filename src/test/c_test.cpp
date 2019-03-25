#include <Rcpp.h>
using namespace Rcpp;

double sample_sd(NumericVector X) {
  double n = X.size();
  double sd1 = sqrt((n-1)/n)*sd(X);
  return sd1;
}

// [[Rcpp::export]]
double pearsoncoeff(NumericVector X, NumericVector Y)
{
  
  return (mean(X*Y)-mean(X)*mean(Y)) / (sample_sd(X)* sample_sd(Y));
}

// [[Rcpp::export]]
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

// [[Rcpp::export]]
float update_score(NumericMatrix X,LogicalVector left_idx_bool,LogicalVector right_idx_bool) {
  float score_nosplit = 0;
  float score_left = 0;
  float score_right = 0;
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
  // Rprintf("%f %f %f\n",score_left,score_right,score_nosplit);
  return score_left+score_right-score_nosplit;
}


// [[Rcpp::export]]
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

// [[Rcpp::export]]
LogicalVector test1(NumericMatrix X, NumericMatrix Y,IntegerVector groups, LogicalVector feature_remaining) {
  int a = feature_remaining.size();
  IntegerVector best_split(X.nrow());
  IntegerVector groups_unique = unique(groups);
  groups_unique = groups_unique[groups_unique!=-1];
  int current_branch = 0;
  LogicalVector current_split_group = groups== current_branch;
  int feature_idx=0;
  NumericVector current_feature_vals = X(_,feature_idx);
  NumericVector valid_vals = current_feature_vals[current_split_group];  
  double current_val = valid_vals[0];
  LogicalVector left_idx_bool = (current_feature_vals <= current_val) & current_split_group;

  LogicalVector right_idx_bool = (current_feature_vals>current_val) & current_split_group;
  float current_split_score = update_score(Y, left_idx_bool, right_idx_bool);
  Rprintf("%f\n",current_split_score);
  return left_idx_bool;

}