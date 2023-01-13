# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 00:45:58 2022

@author: xinyu
"""

import numpy as np
from matplotlib import pyplot as plt
import pickle
from scipy.linalg import sqrtm
from scipy.io import loadmat

dataname = 'TFIMqutip.pkl'
f = open(dataname,'rb')
h,gamma,spin_num,Rs,Ents,tlist, Mats,zs,xs = pickle.load(f)
f.close()

a = np.linspace(0,4990,500,dtype = int)

data1 = loadmat('N4L2R1Jz1h0.mat')
t1 = data1['tlist']
t1 = t1[0][a]
mz1 = data1['mz']
mz1 = mz1[0][a]
mx1 = data1['mx']
mx1 = mx1[0][a]

data4 = loadmat('N4L2R6Jz1h0.mat')
t4 = data4['tlist']
t4 = t4[0][a]
mz4 = data4['mz']
mz4 = mz4[0][a]
mx4 = data4['mx']
mx4 = mx4[0][a]

data8 = loadmat('N4L2R16Jz1h0.mat')
t8 = data8['tlist']
t8 = t8[0][a]
mz8 = data8['mz']
mz8 = mz8[0][a]
mx8 = data8['mx']
mx8 = mx8[0][a]

font1 = {
    'family':'Times New Roman',
    'weight': 'normal',
    'size': 24,
    }

plt.figure(figsize = (12.0,9.0))
plt.tick_params(labelsize=18)         
plt.plot(tlist, zs, '-', color = 'black', label = 'Exact',linewidth=6)
plt.plot(t1,mz1, '-', color = 'blue',label = 'Ansatz II Rank=1',linewidth=4)
plt.plot(t4,mz4, '-', color = 'green',label = 'Ansatz II Rank=6',linewidth=4)
plt.plot(t8,mz8, '-', color = 'red',label = 'Ansatz II Rank = 16',linewidth=4)


plt.legend(prop = font1)

plt.title('Average magnetization of XYZ model ',fontsize=24)
plt.xlabel('t', fontsize = 24)
plt.ylabel('Average magnetization', fontsize = 24)
figname = 'mzXYZL2.png'
plt.savefig(figname)



font1 = {
    'family':'Times New Roman',
    'weight': 'normal',
    'size': 22,
    }

plt.figure(figsize = (12.0,9.0))
plt.tick_params(labelsize=18)         
plt.plot(tlist, xs, '-', color = 'black', label = 'Exact',linewidth=6)
plt.plot(t1,mx1, '-', color = 'blue',label = 'Ansatz II Rank=1',linewidth=4)
plt.plot(t4,mx4, '-', color = 'green',label = 'Ansatz II Rank=6',linewidth=4)
plt.plot(t8,mx8, '-', color = 'red',label = 'Ansatz II Rank=16',linewidth=4)
plt.ylim([-0.1,0.1])
plt.legend(prop = font1)

# plt.plot(gamma*tlist, rho2, marker='.')
plt.title('Average <$\sigma^x$> of XYZ model ',fontsize=24)
plt.xlabel('t', fontsize = 24)
plt.ylabel('Average <$\sigma^x$>', fontsize = 24)
figname = 'mxXYZL2.png'
plt.savefig(figname)
