# Copyright © 2021 LVCS. All Rights Reserved
import numpy as np
import pandas as pd

# 参数设定
string_limit = 700
pcc_limit = 0.5

# STRINGdb转换部分
# STRINGdb数据库导入
stringdb = pd.read_csv("link.csv", delimiter=" ")
# 删除前部的9606人类代号
stringdb["protein1"] = stringdb["protein1"].str.slice(start=5)
stringdb["protein2"] = stringdb["protein2"].str.slice(start=5)
# 筛选出高等可信度数据
stringdb = stringdb[stringdb["combined_score"] >= string_limit]
# 重新排序
stringdb = stringdb.reset_index(drop=True)

# 导入ENSP-序号对应文件
dict1 = pd.read_csv("ensp.csv")
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

string_bool = string_bool == 1

# 各阶段DNB数组初始化
dnb = pd.DataFrame()
range_list = [1, 2, 4, 5, 6]

# DNB计算部分
for j in range_list:

    # 数据导入
    frame = pd.read_csv(str(j) + ".csv")
    del frame["symbol"]
    frame = frame.to_numpy()

    # DNB计算
    max_item = []
    pccin = []
    pccout = []
    sdin = []
    dnb_max_item = []

    # 计算各项的标准差
    sd = []
    for i in range(len(frame)):
        sd.append(np.std(frame[i], ddof=1))

    # 根据数据维度大小可以考虑下面这种计算标准差的方法
    # sd = np.std(frame, axis=1, ddof=1)

    # 计算各项的相关系数
    pc = np.corrcoef(frame)
    pc = np.abs(pc)
    # PCC阈值
    pc_bool = pc >= pcc_limit
    # 根据PCC阈值再次筛选
    pc_bool = string_bool * pc_bool
    pc = pc * pc_bool

    # 计算各项PCCin和SDin
    for i in range(len(pc)):
        pccin_item = np.sum(pc[i, :]) / np.sum(string_bool[i, :])
        sdin_item = np.sum(sd * string_bool[i, :]) / np.sum(string_bool[i, :])
        if str(pccin_item) == "nan":
            pccin_item = 0
        if str(sdin_item) == "nan":
            sdin_item = 0
        pccin.append(pccin_item)
        sdin.append(sdin_item)

    pccin = np.array(pccin)
    sdin = np.array(sdin)

    # 计算各项PCCout
    for i in range(len(pc)):
        x = np.where(pc[i] != 0)[0]
        count = 0
        pccout_item = 0
        for ii in x:
            y = pc[ii] * (string_bool[i] == 0)
            count += np.sum(y != 0)
            pccout_item += np.sum(y)
        if count == 0:
            pccout.append(0)
        else:
            pccout.append(pccout_item / count)

    # 对各阶段数据进行拼接
    dnb_pros = pd.DataFrame(pccin * sdin / pccout)
    dnb_pros.columns = [j]
    dnb = pd.concat([dnb, dnb_pros], axis=1)

symbol = pd.DataFrame(pd.read_csv("1.csv")["symbol"])
dnb = pd.concat([symbol, dnb], axis=1)
dnb = dnb.fillna(0)
dnb.to_csv(
    "DNB-STRING-" + str(string_limit) + "-PCC-" + str(pcc_limit * 1000) + ".csv",
    index=False,
)
