# Copyright © 2021 LVCS. All Rights Reserved
import heapq
import sys

import numpy as np
import pandas as pd

# 文件号
num = str(sys.argv[1])
# 计数阈值
count_num = int(sys.argv[2])
# 输出前N个DNB最大基因
output_num = int(sys.argv[3])

# 初始化
max_item = []
pccin = []
pccout = []
sdin = []
dnb_max_item = []

# 数据导入
frame = pd.read_csv(num + ".csv")
# 构造symbol字典
dict1 = dict(zip(frame.reset_index().values[:, 0], frame.values[:, 0]))

# 删除多余数据
del frame["symbol"]
frame = frame.to_numpy()

# 计算总表格皮尔森相关系数
pc = np.corrcoef(frame)
pc = np.abs(pc)

# 删除对角线
for i in range(len(pc)):
    pc[i, i] = 0

# 找到每项前*个最大相关项
for i, item in enumerate(pc):
    max_item.append(heapq.nlargest(count_num, range(len(item)), item.take))
    pccin_item = np.sum(item)
    sdin_item = np.std(item, ddof=1)
    pccin.append(pccin_item)
    sdin.append(sdin_item)

pccin = np.array(pccin)
sdin = np.array(sdin)

for i in range(len(pc)):
    pccout_item = (np.sum(pccin[max_item[i]]) - pccin[i]) / (count_num - 1)
    pccout.append(pccout_item)

pccout = np.array(pccout)
dnb = pccin * sdin / pccout
dnb_max_item = dnb.argsort()[::-1][0:output_num]
dnb_max_item = dnb_max_item.tolist()

# 输出
print(num, ":")
print("symbol:")
for i in range(output_num):
    print(dict1[dnb_max_item[i]])
print("ID:")
for i in range(output_num):
    print(dnb[dnb_max_item[i]])
