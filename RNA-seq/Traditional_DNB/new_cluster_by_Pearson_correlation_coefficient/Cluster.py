# Copyright © 2021 LVCS. All Rights Reserved
import copy
import sys

import numpy as np
import pandas as pd
from scipy.special import comb

# 文件号
num = str(sys.argv[1])

# 阈值
pccin_min = float(sys.argv[2])
element_limit = int(sys.argv[3])

# 预设聚类参数
cluster_num = int(sys.argv[4])
cluster_len = int(sys.argv[5])

# 数据导入
frame = pd.read_csv(num + ".csv")

# 删除多余数据
del frame["symbol"]
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
# 实例为cluster_num组以内聚类*cluster_len个以内聚类内元素
cluster = np.zeros([cluster_num, cluster_len], dtype=int)

# 计数
i = 0
count = len(pc) - 1
end_flag = 0

while 1:
    # 标记矩阵中无剩余元素，循环跳出，将所有聚类结果输出
    if end_flag == 1:
        # 删除多余全零行列
        cluster = cluster[~(cluster == 0).all(1)]
        idx = np.argwhere(np.all(cluster[..., :] == 0, axis=0))
        cluster = np.delete(cluster, idx, axis=1)
        # 聚类结果输出
        np.savetxt("result" + num + ".csv", cluster, delimiter=",", fmt="%s")
        print("done!")
        break

    # 求出最大相关系数，最大值为零则跳出
    x = np.max(pc_find)

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
    count = count - 2
    pccin = pc[max_peer[0][0], max_peer[1][0]]

    while 1:
        # 若无剩余元素，则跳出
        if count == 0:
            end_flag = 1
            break

        # 极个别情况，仅剩一个元素的处理
        if count == 1:
            if times >= element_limit:
                cluster[i + 1, 0] = np.setdiff1d(range(1,
                                                       len(pc) - 1),
                                                 cluster)[0]
            else:
                cluster[i, times] = np.setdiff1d(range(1,
                                                       len(pc) - 1),
                                                 cluster)[0]
            end_flag = 1
            break

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

        # 极个别情况，仅剩一个元素的处理
        if count == 1:
            if times >= element_limit:
                cluster[i + 1, 0] = np.setdiff1d(range(1,
                                                       len(pc) - 1),
                                                 cluster)[0]
            else:
                cluster[i, times] = np.setdiff1d(range(1,
                                                       len(pc) - 1),
                                                 cluster)[0]
            end_flag = 1
            break

        # 聚类内阈值检测，跳出
        if pccin_ave < pccin_min and times >= element_limit:
            break
        if times > cluster_len - 1:
            break

        # 检测剩余元素，剩余两个直接写入并跳出
        # if len(np.where(pc_find != 0)) == 2:
        if count == 2:
            x = np.max(pc_find)
            last_two = np.where(pc_find == x)
            cluster[i, times] = last_two[0][0]
            cluster[i, times + 1] = last_two[1][0]
            end_flag = 1
            break

    i = i + 1
