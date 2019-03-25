import numpy as np
import pandas as pd
from sklearn.tree import DecisionTreeRegressor

init_group_num = 8
max_partition_level = 3
cor_cutoff = 0.9
min_divide_size = 3
min_group_size = 2
max_iter = 5

def build_regression_tree(X,y):
  regr_1 = DecisionTreeRegressor(max_depth=max_partition_level-1,min_samples_split = min_divide_size)
  regr_1.fit(X, y)
  return regr_1

def assign_tf(tf_data,gene_data,gene_group_table):
    groups_list = []
    group_labels = range(int(max(gene_group_table))+1)
    tf_group_table = np.zeros((len(group_labels),gene_data.shape[1]))-1
    for group_idx in group_labels:
        gene_idx = gene_group_table == group_idx
        if np.sum(gene_idx) >= min_group_size:
            X = tf_data.transpose()
            y = gene_data[gene_idx,].transpose()
            
            regr_1 = build_regression_tree(X,y)
            tf_group_table[group_idx,] = regr_1.apply(tf_data.transpose())
            groups_list.append(regr_1)
    tf_group_table1 = tf_group_table[~np.any(tf_group_table == -1, axis=1)]
    return groups_list,tf_group_table1

def extract_tree_table(estimator,tf_data,tf_list):
    n_nodes = estimator.tree_.node_count
    children_left = estimator.tree_.children_left
    children_right = estimator.tree_.children_right
    feature = estimator.tree_.feature
    threshold = np.round(estimator.tree_.threshold,4)
    node_depth = np.zeros(shape=n_nodes, dtype=np.int64)
    is_leaves = np.zeros(shape=n_nodes, dtype=bool)
    stack = [(0, -1)]  # seed is the root node id and its parent depth
    node_scope = np.ones((n_nodes,tf_data.shape[1]))
    test_id_list = []
    while len(stack) > 0:
        node_id, parent_depth = stack.pop()
        node_depth[node_id] = parent_depth + 1
        # If we have a test node
        if (children_left[node_id] != children_right[node_id]):
            test_id_list.append(node_id)
            stack.append((children_left[node_id], parent_depth + 1))
            stack.append((children_right[node_id], parent_depth + 1))
            left_scope = tf_data[feature[node_id],:] <= threshold[node_id]
            right_scope = tf_data[feature[node_id],:] > threshold[node_id]
            node_scope[children_left[node_id],:] = np.logical_and(left_scope,node_scope[node_id,:])
            node_scope[children_right[node_id],:] = np.logical_and(right_scope,node_scope[node_id,:])
        else:
            is_leaves[node_id] = True  
    tree_table = pd.DataFrame(columns=range(tf_data.shape[1]))       
    for i in range(n_nodes):
        if not is_leaves[i]:
            test_results = 1-(tf_data[feature[i],:] <= threshold[i])
            test_results[node_scope[i,:]==0] = -1
            tree_table.loc[tf_list[feature[i]]] = test_results
    return tree_table   

def calc_mse(x,labels):
    sse_all = 0
    count = 0
    labels_uniq = set(labels)
    for label in labels_uniq:
        label_idx = labels == label
        if np.sum(label_idx) > 1:
            x1 = x[:,label_idx]
            sse_all += np.sum((x1-x1.mean(axis=1,keepdims=True))**2)
            count += x1.shape[0]*x1.shape[1]
    return sse_all/count