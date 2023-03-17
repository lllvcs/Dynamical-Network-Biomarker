import sys

import numpy as np
import pandas as pd

# 读取string阈值
string_limit = int(sys.argv[1])

# stringlink文件读取
stringdb = pd.read_csv("../data/link.csv", delimiter=" ")
stringdb["protein1"] = stringdb["protein1"].str.slice(start=5)
stringdb["protein2"] = stringdb["protein2"].str.slice(start=5)
stringdb = stringdb[stringdb["combined_score"] >= string_limit]
stringdb = stringdb.reset_index(drop=True)

# 构建字典
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

# 文件名列表
namelist = [
    "h000",
    "h005",
    "h012",
    "h021",
    "h029",
    "h036",
    "h045",
    "h053",
    "h060",
    "h069",
    "h077",
    "h084",
    "h093",
    "h101",
    "h108",
]


# 熵计算
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


# 原始样本/baseline计算
origin_frame = pd.read_csv("../data/hbaseline.csv")
del origin_frame["symbol"]
origin_frame = origin_frame.to_numpy()
origin_sd = np.std(origin_frame, ddof=1, axis=1)
origin_pc = np.corrcoef(origin_frame)
origin_pc = np.abs(origin_pc)
origin_pc = origin_pc * string_bool

# 原始样本/baseline边计算
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

# 多样本计算
landscape = pd.DataFrame()
for k in namelist:
    print(k, end="\n")
    origin_frame = pd.read_csv("../data/hbaseline.csv")
    del origin_frame["symbol"]
    append_frame = pd.read_csv("../data/" + str(k) + ".csv")
    del append_frame["symbol"]
    append_frame = pd.concat([origin_frame, append_frame], axis=1)
    origin_frame = origin_frame.to_numpy()
    append_frame = append_frame.to_numpy()

    append_pc = np.corrcoef(append_frame)
    append_pc = np.abs(append_pc)
    append_pc = append_pc * string_bool
    append_sd = np.std(append_frame, ddof=1, axis=1)

    edge_append_entropy = []
    edge_append_sd = []
    for i in range(num):
        for j in range(i - 1):
            if string_bool[i, j] != 0:
                edge_append_sd.append(np.abs(
                    (append_sd[i] + append_sd[j]) / 2))
                edge_append_entropy.append(
                    edge_structure_entropy(append_frame, append_pc, i, j))

    edge_append_entropy = np.array(edge_append_entropy)
    edge_append_entropy = np.abs(edge_append_entropy - edge_origin_entropy)
    edge_append_sd = np.array(edge_append_sd)
    edge_append_sd = np.abs(edge_append_sd - edge_origin_sd)

    landscape_pros = pd.DataFrame(edge_append_sd * edge_append_entropy)
    landscape_pros.columns = [k]
    landscape = pd.concat([landscape, landscape_pros], axis=1)

landscape = pd.concat(
    [pd.DataFrame(edge_origin_list, columns=["node1", "node2"]), landscape],
    axis=1)
landscape = landscape.fillna(0)
landscape.to_csv("edge_structure_entropy_logM-1_" + str(string_limit) + ".csv",
                 index=False)
