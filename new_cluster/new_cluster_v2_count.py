# Copyright © 2021 LVCS. All Rights Reserved
import numpy as np
import pandas as pd
import copy

# 阈值
element_limit = 250

# 数据导入
frame = pd.read_csv('in.csv')

# 删除多余数据
del frame['symbol']
frame = frame.to_numpy()

# 计算总表格皮尔森相关系数
pc = np.corrcoef(frame)
pc = np.abs(pc)

# 删除对角线
for i in range(len(pc)):
    pc[i, i] = 0

# 增加全零行列
pc = np.row_stack((np.zeros(len(pc)), pc))
pc = np.column_stack((np.zeros(len(pc)), pc))

# 生成标记矩阵，便于下一最大相关系数的查找
# 深度拷贝，防止影响原矩阵
pc_temp = copy.deepcopy(pc)
pc_find = copy.deepcopy(pc)
pc_find = np.tril(pc_find)

# 聚类分组初始化
cluster = np.zeros([10000, element_limit], dtype=int)

# 计数
i = 0
count = len(pc) - 1

for i in range(10000):

    # 求出最大相关系数
    x = np.max(pc_find)

    # 标记矩阵中无剩余元素，循环跳出，将所有聚类结果输出
    if x == 0:
        cluster = cluster[~(cluster == 0).all(1)]
        idx = np.argwhere(np.all(cluster[..., :] == 0, axis=0))
        cluster = np.delete(cluster, idx, axis=1)

        np.savetxt('out.csv', cluster, delimiter=',', fmt='%d')
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

    # 聚类内元素计数，初始为2
    times = 2

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

        # 聚类内元素计数加一
        times = times + 1
        count = count - 1

        # 检测剩余元素
        # if len(np.where(pc_find != 0)) == 2:
        if count == 2:
            last_two = np.where(pc_find == x)
            cluster[i, times] = last_two[0][0]
            cluster[i, times + 1] = last_two[1][0]
            break

        # 检测剩余元素，剩余两个直接写入并跳出
        if times >= element_limit:
            break
