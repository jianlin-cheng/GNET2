import math
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from kneed import KneeLocator


# x is n by p matrix, labels is a vector of length n
# for each column of x, calculate sum of log likelihood(probability density of normal distribution) based on each label
# return the sum of log likelihood of all columns
def calc_likelihood_score(x,labels):
    x = x.transpose()
    score_total = 0
    labels_uniq = set(labels)
    for label in labels_uniq:
        label_idx = labels == label
        if np.sum(label_idx) > 1:
            x1 = x[:,label_idx]
            rowmean_x1 = np.mean(x1,axis=1,keepdims=True)
            rowvar_x1 = np.var(x1,axis=1,keepdims=True,ddof=1)
            score_total += np.sum(-((x1-rowmean_x1)**2)/(2*rowvar_x1)-0.5*np.log(2*math.pi*rowvar_x1))
    return score_total

# x is a vector of length n, labels is a vector of length n
# calculate sum of log likelihood(probability density of normal distribution) based on each label
def calc_likelihood_vec(x,labels):
    score_total = 0
    labels_uniq = set(labels)
    for label in labels_uniq:
        label_idx = labels == label
        if np.sum(label_idx) > 1:
            x_i = x[label_idx]
            rowmean_x = np.mean(x_i)
            rowvar_x = np.var(x_i,ddof=1)
            score_total += np.sum(-((x_i-rowmean_x)**2)/(2*rowvar_x)-0.5*np.log(2*math.pi*rowvar_x))
    return score_total

# x is a n by p matrox
# return is the mean of absolute correlation between each pair of columns
def calc_correlation(x):
    if x.shape[0]==1 or x.shape[1]==1:
        return 0
    else:
        x_cor = np.corrcoef(x,rowvar=False)
        cor_mean =  np.mean(np.absolute(x_cor[np.triu_indices(x_cor.shape[0],1)]))
        return cor_mean

# x is n by p matrix, labels is a vector of length p
# for each label, calculate the mean of (absolute correlation between each pair of columns)
# return is a vector pf length p, with values as mean from last step of the label in same position
def get_correlation_list(x,labels):
    cor_list = np.zeros(len(labels))
    for i in set(labels):
        sample_idx = labels==i
        cor_list[sample_idx] = calc_correlation(x[:,sample_idx])
    return cor_list

# transform n by p group table to a vector of length p, which is the label of leaves group for each sample
def get_leaf_labels(group_table):
    leaf_label = np.zeros(group_table.shape[1])
    next_label = 0
    for i in range(group_table.shape[0]):
        current_label = group_table.iloc[i:]
        for j in range(2):
            leaf_label[(current_label==j).values[0]] = next_label
            next_label += 1
    return leaf_label

# X is n by p1 matrix, y is n by p2 matrix
# build regression tree based on each split with highest sum of likelihood score
# keep spliting until all features have been used, or maximum level reached, or
# all subgroups have samples le than min_divide_size, or within group correlation ge than cor_cutoff
# return is a n_feature by n matrix, each row is the split determined by the feature in row name
# samples not determined by the feature in same row are masked as -1
def build_regression_tree_baysian(X,y,max_partition_level = 3,cor_cutoff = 0.9,min_divide_size = 3):
    feature_remaining = list(range(X.shape[1]))
    feature_num = 0
    subgroup_indicator = np.zeros(X.shape[0])
    groups = set([0])
    group_table = pd.DataFrame(columns=range(X.shape[0]),dtype=np.int32)
    while feature_num < max_partition_level and len(groups)>0 and len(feature_remaining)>0:
        best_score = -10**6
        best_feature = None
        for i in groups:
            current_split_group = subgroup_indicator == i
            y_current_split = y[current_split_group,:]
            for j in feature_remaining:
                feature_vals = X[:,j][current_split_group]
                divide_vals = feature_vals[feature_vals != np.max(feature_vals)]
                score_nosplit = calc_likelihood_score(y_current_split,np.ones_like(feature_vals))
                for k in divide_vals:
                    subgroup_divide = 1-(feature_vals <= k)
                    score = calc_likelihood_score(y_current_split,subgroup_divide)-score_nosplit
                    if score > best_score:
                        best_score = score
                        best_feature = j
                        subgroup_labels = np.zeros(len(subgroup_indicator))-1
                        subgroup_labels[current_split_group] = subgroup_divide
        group_table.loc[str(best_feature)] = subgroup_labels
        # subgroups that have samples le than min_divide_size, or correlation ge than cor_cutoff would be masked by -1
        subgroup_indicator_new = get_leaf_labels(group_table)
        subgroup_indicator_new[subgroup_indicator==-1] = -1
        subgroup_indicator = subgroup_indicator_new
        for i in range(2):
            new_divide_i = subgroup_labels==i
            if np.sum(new_divide_i) <= min_divide_size:
                subgroup_indicator[new_divide_i] = -1
            elif calc_correlation(y[new_divide_i,:].transpose()) >= cor_cutoff:
                subgroup_indicator[new_divide_i] = -1
        groups = set(subgroup_indicator[subgroup_indicator!=-1])
        feature_remaining = [i for i in feature_remaining if i != best_feature]
        feature_num+=1
    return group_table
    

# Assign TF to groups with baysian
def assign_tf_baysian(tf_data,gene_data,gene_group_table,max_partition_level = 3,
             cor_cutoff = 0.9,min_divide_size = 3,min_group_size = 2):
    group_labels = range(int(max(gene_group_table))+1)
    tf_group_table = np.zeros((len(group_labels),gene_data.shape[1]))-1
    groups_list = []
    for group_idx in group_labels:
        gene_idx = gene_group_table == group_idx
        if np.sum(gene_idx) >= min_group_size:
            X = tf_data.transpose()
            y = gene_data[gene_idx,:].transpose()
            group_table = build_regression_tree_baysian(X,y,max_partition_level,cor_cutoff,min_divide_size)
            groups_list.append(group_table)
            tf_group_table[group_idx,] = get_leaf_labels(group_table)
    # remove groups that with number of genes less than min_group_size
    tf_group_table1 = tf_group_table[~np.any(tf_group_table == -1, axis=1)]
    return groups_list,tf_group_table1
        
# Assign genes to groups 
def assign_gene(gene_data,tf_group_table):
    gene_group_table = np.zeros(gene_data.shape[0])-1
    for gene_idx in range(gene_data.shape[0]):
        exp_gene = gene_data[gene_idx,:]
        def calc_likelihood_vec1(leaf_label):
            return calc_likelihood_vec(exp_gene,leaf_label)
        gene_group_table[gene_idx] = np.argmax(np.apply_along_axis(calc_likelihood_vec1,1,tf_group_table))
    return gene_group_table


def run_gnet(input_data,tf_list,init_group_num = 8,max_partition_level = 3,cor_cutoff = 0.9,min_divide_size = 3,
             min_group_size = 2,max_iter = 5,min_group_num = 3,max_cluster_size = 20):
    gene_data = input_data.drop(tf_list,axis=0).values
    tf_data = input_data.loc[tf_list].values
    gene_names = input_data.drop(tf_list,axis=0).index
    # initial clustering
    print('Determine initial group number')
    avg_cor_list = np.zeros(min(max_cluster_size,gene_data.shape[0])-1)
    for i in range(min(max_cluster_size,gene_data.shape[0])-1):
        gene_group_table = KMeans(n_clusters=i+2, random_state=0).fit(gene_data).labels_
        avg_cor_list[i] = np.mean(get_correlation_list(gene_data.transpose(),gene_group_table))

    # Determine knee point
    if len(avg_cor_list) > min_group_num:
        order = np.argsort(-avg_cor_list)
        y = avg_cor_list[order]
        x = np.array(range(1, len(y)+1))
        kn = KneeLocator(x, y, curve='convex', direction='decreasing')
        if kn.knee is None:
            groups_keep = order
        else:
            groups_keep = order[range(kn.knee)]
    else:
        groups_keep = np.array(range(min_group_num))
    gene_group_table[np.logical_not(np.isin(gene_group_table , groups_keep))] = -1

    # The E-M process
    for i in range(max_iter):
        print('iteration '+str(i))
        groups_list,tf_group_table = assign_tf_baysian(tf_data,gene_data,gene_group_table,max_partition_level,
             cor_cutoff,min_divide_size,min_group_size)
        gene_group_table_new = assign_gene(gene_data,tf_group_table)
        if np.array_equal(gene_group_table_new, gene_group_table):
            exit
        else:
           gene_group_table =  gene_group_table_new

    # formatting output
    group_idx = 0
    tree_table_all = pd.DataFrame(columns=['group','index']+list(range(groups_list[0].shape[1])))
    gene_list_all = pd.DataFrame(columns=['group','gene'])
    for i in range(len(groups_list)):
        tree_table = groups_list[i]
        gene_list = np.array(gene_names)[gene_group_table==i].tolist()
        tree_table_ordered = tree_table.reset_index()
        new_tree_table = pd.concat([pd.DataFrame(np.zeros((tree_table_ordered.shape[0],1))+group_idx,columns=['group']),
                               tree_table_ordered],axis=1, sort=False)
        tree_table_all = pd.concat([tree_table_all,new_tree_table],axis=0,sort=False, ignore_index=True)
        
        new_gene_table = pd.concat([pd.DataFrame(np.zeros((len(gene_list),1))+group_idx,columns=['group']),
                                    pd.DataFrame(gene_list,columns=['gene'])],axis=1, sort=False)
        
        gene_list_all = pd.concat([gene_list_all,new_gene_table],axis=0,sort=False, ignore_index=True)
        group_idx += 1
    tree_table_all.columns = ['group','feature']+ input_data.columns.values.tolist()

    return gene_list_all,tree_table_all
 
    
#test
#rnaseq_file = '/home/chen/Dropbox/MU/workspace/new/data/test/rnaseq_data.csv'
#tf_list_file  = '/home/chen/Dropbox/MU/workspace/new/data/test/tf_list.csv'

rnaseq_file = 'C:/Users/chen/Documents/Dropbox/MU/workspace/new/data/test/rnaseq_data.csv'
tf_list_file  = 'C:/Users/chen/Documents/Dropbox/MU/workspace/new/data/test/tf_list.csv'

init_group_num = 12
max_partition_level = 3
cor_cutoff = 0.9
min_divide_size = 3
min_group_size = 2
max_iter = 5
min_group_num = 3
max_cluster_size = 20

input_data = pd.read_csv(rnaseq_file,index_col=0)
with open(tf_list_file) as f:
    tf_list = f.read().splitlines()
        
gene_list_all,tree_table_all = run_gnet(input_data,tf_list,init_group_num,max_partition_level,cor_cutoff,
                                        min_divide_size,min_group_size,max_iter,min_group_num,max_cluster_size)  

#input_data = pd.read_csv(rnaseq_file,index_col=0)
#with open(tf_list_file) as f:
#    tf_list = f.read().splitlines()
#
#gene_data = input_data.drop(tf_list,axis=0).values
#tf_data = input_data.loc[tf_list].values
#    
#gene_group_table=np.concatenate([np.zeros(1500),np.zeros(2000)+1,np.zeros(2000)+2,np.zeros(1000)+3])
#groups_list,tf_group_table = assign_tf_baysian(tf_data,gene_data,gene_group_table)
#avg_cor_list = np.zeros(tf_group_table.shape[0])
#for i in range(tf_group_table.shape[0]):
#    x_group = gene_data[gene_group_table==i,]
#    avg_cor_list[i] = np.mean(get_correlation_list(x_group,tf_group_table[i,:]))
#groups_keep = np.argsort(-avg_cor_list)

