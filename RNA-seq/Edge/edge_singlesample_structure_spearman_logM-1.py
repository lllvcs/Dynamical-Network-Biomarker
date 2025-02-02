import sys

import numpy as np
import pandas as pd
import scipy.stats as stats

string_limit = int(sys.argv[1])

stringdb = pd.read_csv("../data/link.csv", delimiter=" ")
stringdb["protein1"] = stringdb["protein1"].str.slice(start=5)
stringdb["protein2"] = stringdb["protein2"].str.slice(start=5)
stringdb = stringdb[stringdb["combined_score"] >= string_limit]
stringdb = stringdb.reset_index(drop=True)

dict1 = pd.read_csv("../data/ensp.csv")
length = len(np.unique(dict1["symbol"].values))

d = {dict1["ensp"][i]: int(dict1["id"][i]) for i in range(len(dict1))}
d2 = {int(dict1["id"][i]): dict1["symbol"][i] for i in range(len(dict1))}
string_bool = np.zeros([length, length])

for i in range(len(stringdb)):
    if (stringdb["protein1"][i] in d) & (stringdb["protein2"][i] in d):
        string_bool[d[stringdb["protein1"][i]], d[stringdb["protein2"][i]]] = 1

np.fill_diagonal(string_bool, 0)
string_bool = string_bool == 1


def edge_structure_entropy(data, pc, p1, p2):
    left_not_zero = np.setdiff1d(np.where(pc[p1] != 0), p2)
    right_not_zero = np.setdiff1d(np.where(pc[p2] != 0), p1)
    len_left_not_zero = len(left_not_zero)
    len_right_not_zero = len(right_not_zero)

    if len_left_not_zero < 2 and len_right_not_zero < 2:
        return 0
    if len_left_not_zero < 2 and len_right_not_zero > 1:
        right_prob = (abs(pc[p2][right_not_zero])) / np.sum(
            abs(pc[p2][right_not_zero]))
        right_entropy = -np.sum(
            np.average(data[right_not_zero, :], axis=1) * right_prob *
            np.log2(right_prob)) / np.log2(len_right_not_zero)
        entropy = right_entropy / 2
        return entropy
    if len_left_not_zero > 1 and len_right_not_zero < 2:
        left_prob = (abs(pc[p1][left_not_zero])) / np.sum(
            abs(pc[p1][left_not_zero]))
        left_entropy = -np.sum(
            np.average(data[left_not_zero, :], axis=1) * left_prob *
            np.log2(left_prob)) / np.log2(len_left_not_zero)
        entropy = left_entropy / 2
        return entropy

    left_prob = (abs(pc[p1][left_not_zero])) / np.sum(
        abs(pc[p1][left_not_zero]))
    right_prob = (abs(pc[p2][right_not_zero])) / np.sum(
        abs(pc[p2][right_not_zero]))
    left_entropy = -np.sum(
        np.average(data[left_not_zero, :], axis=1) * left_prob *
        np.log2(left_prob)) / np.log2(len_left_not_zero)
    right_entropy = -np.sum(
        np.average(data[right_not_zero, :], axis=1) * right_prob *
        np.log2(right_prob)) / np.log2(len_right_not_zero)
    entropy = (left_entropy + right_entropy) / 2
    return entropy


origin_frame = pd.read_csv("../data/hbaseline.csv")
del origin_frame["symbol"]
origin_frame = origin_frame.to_numpy()
origin_sd = np.std(origin_frame, ddof=1, axis=1)
origin_pc = stats.spearmanr(origin_frame, axis=1).correlation
origin_pc = np.abs(origin_pc)
origin_pc = origin_pc * string_bool

edge_origin_list = []
edge_origin_entropy = []
edge_origin_sd = []

num = len(origin_frame)
for i in range(num):
    for j in range(i - 1):
        if string_bool[i, j] != 0:
            edge_origin_list.append([d2[i], d2[j]])
            edge_origin_entropy.append(
                edge_structure_entropy(origin_frame, origin_pc, i, j))
            edge_origin_sd.append((origin_sd[i] + origin_sd[j]) / 2)
edge_origin_list = np.array(edge_origin_list)
edge_origin_entropy = np.array(edge_origin_entropy)
edge_origin_sd = np.array(edge_origin_sd)

landscape = pd.DataFrame()

single_sample = pd.read_csv("../data/exprSet_new.csv", index_col=0)
num_sample = single_sample.columns
num_gene = len(single_sample.iloc[:, 0])

for k, item in enumerate(num_sample):
    print(k, end="\n")
    new = np.hstack((origin_frame,
                     np.reshape(single_sample.iloc[:, k].values,
                                (num_gene, 1))))

    new_pc = stats.spearmanr(new, axis=1).correlation
    new_pc = np.abs(new_pc)
    new_pc = new_pc * string_bool
    new_sd = np.std(new, ddof=1, axis=1)

    edge_append_entropy = []
    edge_append_sd = []

    for i in range(num):
        for j in range(i - 1):
            if string_bool[i, j] != 0:
                edge_append_sd.append(np.abs((new_sd[i] + new_sd[j]) / 2))
                edge_append_entropy.append(
                    edge_structure_entropy(new, new_pc, i, j))

    edge_append_entropy = np.array(edge_append_entropy)
    edge_append_entropy = np.abs(edge_append_entropy - edge_origin_entropy)
    edge_append_sd = np.array(edge_append_sd)
    edge_append_sd = np.abs(edge_append_sd - edge_origin_sd)

    landscape_pros = pd.DataFrame(edge_append_sd * edge_append_entropy)
    landscape_pros.columns = [item]
    landscape = pd.concat([landscape, landscape_pros], axis=1)

landscape = pd.concat(
    [pd.DataFrame(edge_origin_list, columns=["node1", "node2"]), landscape],
    axis=1)
landscape = landscape.fillna(0)
landscape.to_csv(
    "edge_singlesample_structure_entropy_spearman_logM-1_" +
    str(string_limit) + ".csv",
    index=False,
)
