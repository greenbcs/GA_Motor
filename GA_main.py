# -*- coding:utf-8 -*-
import time
start = time.clock()
import numpy as np
import nl_GA_funs as fun
reload(fun)
import matplotlib.pyplot as plt
####所有参数
popsize = 300 #种群数量
pc = 0.7  #交叉概率
pm = 0.002 #变异概率
runtime = 300 #代码运行次数
ddl = int(0.45*runtime) #用于控制退出时间。
results = [[]]#用于存放结果
bestindividual = [] #用于存放最好的结果
bestfit = 0 #用来存放最好的拟合度
fitvalue = [] #用来存放所有的拟合程度
tempop = [[]] #用于存放历史种群?
#bound = [[0,0.9*np.pi],[0,0.9*np.pi],[0,0.9*np.pi],[0,0.9*np.pi],[0,0.9*np.pi]]#用来存放约束，#表示函数1的每个自变量的范围
bound = [[-1000,1000],[-1000,1000]]#用来存放约束，#表示函数1的每个自变量的范围
chromlength = len(bound) #自变量的数目,染色体长度，函数1的染色体长度为5，函数2为2
######以下内容不需要修改######
forpop = np.zeros([popsize,chromlength])
for lie in range(len(bound)):
	forpop[:,lie] = np.random.uniform(bound[lie][0],bound[lie][1],(1,popsize)) #用来初始化种群，1是全部x的下限，2是全部x的上线，括号里面是数量。
pop = forpop.copy()
#到这里我们完成了初始种群的生成
minists = [] #用于记录每一届的最优值
allmin = 100 #用于记录全局最小值
plt.figure()
for t in range(runtime):
	ys = fun.ys(pop) #根据生成的初始种群计算
	fits = 1.0/np.array(ys)  #计算适应度，因为是求极小值，所以取倒数
	minist = min(ys).copy() #进化前的最小值
	#记录此时最优值的位置和结果
	ministindex = ys.index(minist)
	ministpop = pop[ministindex]
	sumf = sum(fits)
	best = max(fits)
	idv_pro = fits/sumf
	#######开始建立新的种群##########
	ret = [] #用于存放新的种群
	for i in range(len(pop)):#这个之所以要循环这么多次，是想要保持，种群大小不变
		pick = np.random.uniform(0,1)
		while pick==0:
			pick = np.random.uniform(0, 1)
		for j in range(len(pop)): #这里再对每一个随机数进行遍历操作，用于寻找较优解
			pick = pick - idv_pro[j]
			if pick<0:
				ret.append(pop[j].tolist()) #这意味着我们会把适应度较高的元素，多复制几个作为新的种群
				break #在完成了添加之后，也就跳出循环了
	ret = np.array(ret)
	pop = ret
	# ###########以上完成了新的种群的建立##########
	pop = fun.jiaocha(popsize, chromlength, pc, bound, pop)
	pop = fun.bianyi(popsize,chromlength,pop,bound,runtime,t,pm)
	################最差值更新为最优值###############
	ys_1 = fun.ys(pop) #根据生成的初始种群计算
	minist_1 = min(ys_1) #进化后的最小值
	maxist_1 = max(ys_1) #进化后的最大值
	maxistindex_1 = ys_1.index(maxist_1)
	ministindex_1 = ys_1.index(minist_1)
	if minist_1>=allmin:#也就是如果这次的值还不如上次的值好的话
		print minist_1>=minist
		pop[ministindex_1] = allminpop.copy() #就把全剧最佳给他
		minist_1=minist
	pop[maxistindex_1] = pop[ministindex_1].copy()
	if t>1:
		plt.plot([t,t-1],[minist_1,minists[-1]])
	minists.append(minist_1)
	if minist_1<allmin:
		ret0 = pop.copy()
		allmin = minist_1
		allminpop = ret0[ys_1.index(minist_1)].copy()
	if t>ddl and minist_1 >= minists[-ddl+1]: #ddl用于控制时间
		break
	print t
print '最优解：'+str(allminpop)
print '最优适应度：'+str(allmin)
plt.show()
end = time.clock()
print '所用时间：'+str(end-start)