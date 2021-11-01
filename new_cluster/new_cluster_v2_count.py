import numpy as np
import pandas as pd
import copy

# 聚类内元素阈值
element_limit = 250

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

# 增加全零行列
pc = np.row_stack((np.zeros(len(pc)), pc))
pc = np.column_stack((np.zeros(len(pc)), pc))

# 生成标记矩阵，便于下一最大相关系数的查找
# 深度拷贝，防止影响原矩阵
pc_temp = copy.deepcopy(pc)
pc_find = copy.deepcopy(pc)
pc_find = np.tril(pc_find)

# 聚类分组初始化
cluster = np.empty([10000, element_limit], dtype=int)

i = 0
for i in range(10000):
    # 求出最大相关系数
    x = np.max(pc_find)

    # 标记矩阵中无剩余元素，循环跳出，将所有聚类结果输出
    if x == 0:
        np.savetxt('result50-1_cluster_250.csv', cluster, delimiter=',')
        break

    # 查找最大相关系数的坐标
    max_peer = np.where(pc_temp == np.max(pc_temp))

    # 最大相关系数对写入聚类分组
    cluster[i][0] = max_peer[0][0]
    cluster[i][1] = max_peer[1][0]

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
        group = 0

        # 聚类内元素对应相关系数行相加
        for j in range(times):  # j
            group = group + pc_temp[cluster[i][j]]
        group = group / times

        # 查找外部相关系数最大的元素
        max_peer_next = np.where(group == np.max(group))[0][0]

        # 将该元素加入聚类内
        cluster[i][times] = max_peer_next

        # 在标记矩阵中移除该加入聚类的元素
        pc_find[max_peer_next] = 0
        pc_find[:, max_peer_next] = 0
        pc_temp[:, max_peer_next] = 0

        # 聚类内元素计数加一
        times = times + 1

        # 聚类内元素达到阈值，跳出
        if times >= element_limit:
            break
