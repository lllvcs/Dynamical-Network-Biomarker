# Copyright © 2021 LVCS. All Rights Reserved
import numpy as np
import pandas as pd
from scipy.special import comb
import copy

# 阈值
pccin_min = 0.75
element_limit = 20

# 数据导入
frame = pd.read_csv('read.csv')

# 删除多余数据
del frame['symbol']
frame = frame.to_numpy()

# 计算总表格皮尔森相关系数
pc = np.corrcoef(frame)
pc = np.abs(pc)

# 删除对角线
for i in range(len(pc)):
    pc[i, i] = 0
pc = np.row_stack((np.zeros(len(pc)), pc))
pc = np.column_stack((np.zeros(len(pc)), pc))

# 深度拷贝，防止影响原矩阵
pc_temp = copy.deepcopy(pc)
pc_find = copy.deepcopy(pc)
pc_find = np.tril(pc_find)

# 聚类分组初始化，根据情况设置预留空间
cluster = np.empty([10000, 1000], dtype=int)
count = len(pc) - 1

# 聚类计数
i = 0

# 求出最大相关系数
x = np.max(pc_find)

for i in range(10000):

    # 求出最大相关系数，最大值为零则跳出
    x = np.max(pc_find)

    # 标记矩阵中无剩余元素，循环跳出，将所有聚类结果输出
    if x == 0:
        cluster = cluster[~(cluster == 0).all(1)]
        idx = np.argwhere(np.all(cluster[..., :] == 0, axis=0))
        cluster = np.delete(cluster, idx, axis=1)

        np.savetxt('result.csv', cluster, delimiter=',', fmt='%d')
        print("done!")
        break

    # 查找最大相关系数的坐标
    max_peer = np.where(pc_find == x)

    # 最大相关系数对写入聚类分组
    cluster[i, 0] = max_peer[0][0]
    cluster[i, 1] = max_peer[1][0]

    # 在标记矩阵中移除写入聚类分组的数据
    pc_find[max_peer[0][0]] = 0
    pc_find[max_peer[1][0]] = 0
    pc_find[:, max_peer[0][0]] = 0
    pc_find[:, max_peer[1][0]] = 0
    pc_temp[:, max_peer[0][0]] = 0
    pc_temp[:, max_peer[1][0]] = 0

    # 聚类内元素计数，初始为2，PCCin即为相关系数
    times = 2
    pccin = pc[max_peer[0][0], max_peer[1][0]]

    while 1:

        # 初始化
        group = 0

        # 聚类内元素对应相关系数行相加
        for j in range(times):  # j
            group = group + pc_temp[cluster[i, j]]
        group = group / times

        # 查找外部相关系数最大的元素
        max_peer_next = np.where(group == np.max(group))[0][0]

        # 将该元素加入聚类内
        cluster[i, times] = max_peer_next

        # 在标记矩阵中移除该加入聚类的元素
        pc_find[max_peer_next] = 0
        pc_find[:, max_peer_next] = 0
        pc_temp[:, max_peer_next] = 0

        # 计算聚类内PCCin之和，判断是否跳出
        for j in range(times):
            pccin = pccin + pc[(cluster[i, j])][(cluster[i, times])]

        # 聚类内元素计数加一
        times = times + 1
        count = count - 1

        # 计算聚类内PCCin均值
        pccin_ave = pccin / comb(times, 2)

        # 聚类内阈值检测，跳出
        if pccin_ave < pccin_min and times >= element_limit:
            break
        if times >= 1000:
            break

        # 检测剩余元素
        # if len(np.where(pc_find != 0)) == 2:
        if count == 2:
            last_two = np.where(pc_find == x)
            cluster[i, times] = last_two[0][0]
            cluster[i, times+1] = last_two[1][0]
            break
