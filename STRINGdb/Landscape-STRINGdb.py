# Copyright © 2021 LVCS. All Rights Reserved
import numpy as np
import pandas as pd

# STRINGdb转换部分
# STRINGdb数据库导入
stringdb = pd.read_csv("link.csv", delimiter=" ")
# 删除前部的9606人类代号
stringdb["protein1"] = stringdb["protein1"].str.slice(start=5)
stringdb["protein2"] = stringdb["protein2"].str.slice(start=5)
# 筛选出高等可信度数据
stringdb = stringdb[stringdb["combined_score"] >= 700]
# 重新排序
stringdb = stringdb.reset_index(drop=True)

# 导入ENSP-序号对应文件
dict1 = pd.read_csv('ensp.csv')
length = len(np.unique(dict1["ENSG"].values))

# 构建字典文件
d = {dict1["ENSP"][i]: int(dict1["id"][i]) for i in range(len(dict1))}

string_bool = np.zeros([length, length])

# STRINGdb转换为邻接矩阵
for i in range(len(stringdb)):
    if (stringdb["protein1"][i] in d) & (stringdb["protein2"][i] in d):
        string_bool[d[stringdb["protein1"][i]], d[stringdb["protein2"][i]]] = 1

for i in range(len(string_bool)):
    string_bool[i, i] = 0

string_bool = (string_bool == 1)

# Landscape计算部分
for j in range(2, 7):
    # 数据导入
    origin_frame = pd.read_csv('1.csv')
    del origin_frame['symbol']
    append_frame = pd.read_csv(j+'.csv')
    del append_frame['symbol']
    append_frame = pd.concat([origin_frame, append_frame], axis=1)
    origin_frame = origin_frame.to_numpy()
    append_frame = append_frame.to_numpy()

    # 标准差计算
    origin_sd = []
    append_sd = []
    for i in range(len(origin_frame)):
        origin_sd.append(np.std(origin_frame[i], ddof=1))
        append_sd.append(np.std(append_frame[i], ddof=1))
    origin_sd = np.array(origin_sd)
    append_sd = np.array(append_sd)
    append_sd = np.abs(append_sd - origin_sd)

    # 相关系数计算
    origin_pc = np.corrcoef(origin_frame)
    origin_pc = np.abs(origin_pc)
    origin_pc = origin_pc*string_bool
    append_pc = np.corrcoef(append_frame)
    append_pc = np.abs(append_pc)
    append_pc = append_pc*string_bool

    # 熵计算
    delta_entropy = []
    for i in range(len(origin_pc)):
        x = np.where(origin_pc[i] != 0)[0]
        y = np.where(append_pc[i] != 0)[0]
        entropy = np.abs(np.sum(np.log2(append_pc[i][y]/np.sum(append_pc[i])))/np.log2(len(y)) \
        - np.sum(np.log2(origin_pc[i][x]/np.sum(origin_pc[i])))/np.log2(len(x)))
        if str(entropy) == 'nan':
            entropy = 0
        delta_entropy.append(entropy)

    # Landscape输出
    landscape = append_sd * delta_entropy
    np.savetxt("landscape"+j+".csv", delta_entropy, delimiter=',', fmt='%s')