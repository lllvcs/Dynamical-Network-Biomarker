import numpy as np
import pandas as pd
from scipy.special import perm, comb
import random
import copy

# 数据导入
frame1 = pd.read_csv('test.csv')
frame2 = copy.deepcopy(frame1)
# 删除多余数据
del frame1['symbol']
frame1 = frame1.to_numpy()

# 计算总表格皮尔森相关系数
pc = np.corrcoef(frame1)
pc = np.abs(pc)
# 删除对角线
for i in range(len(pc)):
    pc[i, i] = 0
pc = np.row_stack((np.zeros(len(pc)), pc))
pc = np.column_stack((np.zeros(len(pc)), pc))
# 深度拷贝，防止影响原矩阵
pc_temp = copy.deepcopy(pc)

# 聚类分组初始化
cluster = np.empty([10000, 1000], dtype=int)

i = 0

for i in range(10000):
    # 求出最大相关系数
    x = np.max(pc_temp)
    print("i:", i)
    if i == 102:
        print("hhh")

    if x == 0:
        break

    max_peer = np.where(pc_temp == np.max(pc_temp))
    cluster[i][0] = max_peer[0][0]
    cluster[i][1] = max_peer[1][0]
    times = 0
    while 1:
        group = 0

        #pccin = 0
        #pcc_ave = 0
        count = np.bincount(cluster[i])[1:].sum()  # 有BUG！会卡在第0位无限循环
        print("count:", count)
        for j in range(count):
            for k in range(count):
                pc_temp[(cluster[i][j])][(cluster[i][k])] = 0
                #pccin = pccin + pc[(cluster[i][j])][(cluster[i][k])]
                # print("pccin:", pccin)
        #pcc_ave = (pccin / 2) / (comb(count, 2))
        #print("pcc_ave:", pcc_ave)

        for j in range(count):  # j
            # print("cluster member:", cluster[i][j])
            # print("i,j:", i, j)
            group = group + pc_temp[cluster[i][j]]

        group = group / count
        print("max_corr", group.max())
        max_peer_next = np.where(group == np.max(group))[0][0]
        print("next peer:", max_peer_next)
        cluster[i][count] = max_peer_next
        times = times + 1
        print("times:", times, end="\n\n")

        if times > 250:
            print("break")
            break
'''
        if pcc_ave < 0.9:
            print("break")
            break
'''




