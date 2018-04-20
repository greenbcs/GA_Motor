# -*- coding:utf-8 -*-
#This file was for save functions to run GA_main
#Author : Zherui Liu
#Date : 2018-3-26
import numpy as np

def ys(pop): #pop is a multi-array\
    ys=[]
    for i in range(len(pop)):
        #函数1
        # y = -5 * np.sin(pop[i][0]) * np.sin(pop[i][1]) * np.sin(pop[i][2]) * np.sin(pop[i][3]) * np.sin(pop[i][4]) - np.sin(5*pop[i][0]*5*pop[i][1]*5*pop[i][2]*5*pop[i][3]*5*pop[i][4]) + 8
        #函数2
        y = -20 * np.exp(-0.2 * ((pop[i][0] ** 2 + pop[i][1] ** 2) / 2) ** (1 / 2.0)) - np.exp(
            1 / 2.0 * (np.cos(2 * np.pi * pop[i][0]) + np.cos(2 * np.pi * pop[i][1]))) + 20 + 2.71289
        ys.append(y)
    return ys

def test(bound,pop,a):
    for n in range(len(pop[a])):
        if pop[a][n] < bound[n][1] and bound[n][0] < pop[a][n]:
            flag1 = 1  # 1表示可行，所以只要有一个不可行为0，那么总的结果就会是不可行
            continue
        else:
            flag1 = 0
            break
    return flag1

def jiaocha(popsize,chromlength,pc,bound,pop):
    ###########以下开始进行交叉操作###############
    for i in range(popsize):
        # print i
        # 随机挑选两个染色体进行交叉操作
        pick = np.random.uniform(0, 1, (1, 2))
        index = np.ceil(popsize * pick) - 1
        # 交叉概率决定是否交叉
        pick = np.random.uniform(0, 1)
        while pick == 0:
            pick = np.random.uniform(0, 1)
        if pick > pc:
            continue  # 这个pick是为了筛选交叉概率,如果>的话，就跳过下面的内容
        # 接下来将要进行的就是开始随机选取交叉位置
        flag = 0
        g = 0
        while flag == 0:
            g = g + 1
            # print g
            pick = np.random.uniform(0, 1)  # 这个pick是用来决定交叉位置
            while pick == 0:
                pick = np.random.uniform(0, 1)
            pos = int(np.ceil(pick * chromlength) - 1)
            # 交叉开始
            pick = np.random.uniform(0, 1)  # 这个pick是用于进行交叉
            a = int(index[0][0])
            b = int(index[0][1])
            v1 = pop[a, pos]
            v2 = pop[b, pos]
            pop[a, pos] = pick * v2 + (1 - pick) * v1
            pop[b, pos] = pick * v1 + (1 - pick) * v2
            # 交叉结束,测试交叉之后的结果是否可行
            for n in range(len(pop[a])):
                if pop[a][n] < bound[n][1] and bound[n][0] < pop[a][n]:
                    flag1 = 1  # 1表示可行，所以只要有一个不可行为0，那么总的结果就会是不可行
                    continue
                else:
                    flag1 = 0
                    break
            flag2 = test(bound, pop, b)
            if flag1 * flag2 == 0:
                flag = 0
            else:
                flag = 1
    return pop

def bianyi(popsize,chromlength,pop,bound,runtime,t,pm):
    ###########变异开始#############
    for i in range(popsize):
        # print i
        # 随机选择一个个个体进行变异
        pick = np.random.uniform(0, 1)
        while pick == 0:
            pick = np.random.uniform(0, 1)
        index = int(pick * popsize)
        # 变异概率决定是否最终变异
        pick = np.random.uniform(0, 1)
        if pick > pm:
            continue
        flag = 0
        while flag == 0:  # 表示如果结果不合理就一直进行这个循环
            # 选取变异位置：
            pick = np.random.uniform(0, 1)
            while pick == 0:
                pick = np.random.uniform(0, 1)
            pos = int(pick * chromlength)
            v = pop[i][pos]
            v1 = v - bound[pos][0]
            v2 = bound[pos][1] - v
            # 开始变异
            pick = np.random.uniform(0, 1)
            if pick > 0.5:
                delta = v2 * (1 - pick ** ((1.0 - t / runtime) ** 2))
                pop[i][pos] = v + delta
            else:
                delta = v1 * (1 - pick ** ((1.0 - t / runtime) ** 2))
                pop[i][pos] = v - delta
            # 变异结束
            flag = test(bound, pop, i)
    return pop

