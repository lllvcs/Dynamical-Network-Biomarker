# -*- coding: UTF-8 -*-
import numpy as np
import pandas as pd
import copy
import random
from scipy.special import comb

dnb = pd.DataFrame()

for iii in range(6):
    print(iii+1)
    # 导入聚类结果的原始数据raw，初始化
    cluster_result_raw = np.loadtxt(
        ("in" + str(iii + 1) + ".csv"), delimiter=",")
    total = np.sum(cluster_result_raw != 0)
    cluster_done = np.zeros(total + 1, dtype=int)
    cluster_len = len(cluster_result_raw)
    
    # 整理聚类
    for i in range(cluster_len):
        for j in range(len(cluster_result_raw[0][:])):
            if cluster_result_raw[i][j] == 0:
                break
            else:
                cluster_done[int(cluster_result_raw[i, j])] = str(i + 1)

    # 删除多余项
    cluster_done = np.delete(cluster_done, cluster_done == 0)

    # 转换为df类型，方便输出
    cluster_done = pd.DataFrame(cluster_done)
    cluster_done.columns = ['cluster']

    # 导入原始基因表达量表格
    origin = pd.read_csv(str(iii + 1) + '.csv')
    origin.insert(1, 'cluster', cluster_done)
    origin = origin.sort_values(by=['cluster', 'symbol'])
    origin.insert(2, 'num', np.linspace(
        1, total, endpoint=True, num=total, dtype=int))
    origin.to_csv('cluster'+str(iii+1)+'.csv',
                  columns=['symbol', 'cluster'], sep=',', index=False, header=True)

    x = copy.deepcopy(origin)
    y = copy.deepcopy(origin)

    # 删除多余数据
    del x['cluster']
    del x['symbol']
    del x['num']
    x = x.values

    # 计算总表格皮尔森相关系数
    pc = np.corrcoef(x)
    pc = np.abs(pc)

    # 各阶段DNB数组初始化
    dnb_pros = np.array([])
    # 计算DNB
    for ii in range(1, int(y.max().cluster + 1)):
        start = y[(y['cluster'] == ii)].min().num
        end = y[(y['cluster'] == ii)].max().num
        pccin = np.triu(pc[start - 1:end, start - 1:end])
        pccin_ave = (np.sum(pccin) - (end - start + 1)) / \
            comb(end - start + 1, 2)
        sdin = 0
        for i in range(start - 1, end):
            sdin = sdin + np.std(x[i, :], ddof=1)
        sdin_ave = sdin / (end - start + 1)
        i = 1
        pccout = 0
        times = 5000
        while i <= times:
            rand = random.randint(0, len(x) - 1)
            if start <= rand <= end:
                continue
            else:
                pccout = pccout + \
                    np.sum(pc[rand, start - 1:end]) / (end - start + 1)
                i = i + 1
        pccout_ave = pccout / times
        dnb_pros = np.append(dnb_pros, pccin_ave * sdin_ave / pccout_ave)

    # 对各阶段数据进行拼接
    dnb_pros = pd.DataFrame(dnb_pros)
    dnb_pros.columns = [iii + 1]
    dnb = pd.concat([dnb, dnb_pros], axis=1)

# 输出DNB总结果
np.savetxt('DNB.csv', dnb, delimiter=',')
