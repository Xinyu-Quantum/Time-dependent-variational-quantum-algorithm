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

data1 = loadmat('N3L3R1Jz1h0.25.mat')

a = np.linspace(0,4990,500,dtype = int)

t1 = data1['tlist']
t1 = t1[0][a]
mz1 = data1['mz']
mz1 = mz1[0][a]
mx1 = data1['mx']
mx1 = mx1[0][a]
data8 = loadmat('N3L3R8Jz1h0.25.mat')
t8 = data8['tlist']
t8 = t8[0][a]
mz8 = data8['mz']
mz8 = mz8[0][a]
mx8 = data8['mx']
mx8 = mx8[0][a]

step1 = np.linspace(0,490,50, dtype = int)

font1 = {
    'family':'Times New Roman',
    'weight': 'normal',
    'size': 24,
    }

plt.figure(figsize = (12.0,9.0))
plt.tick_params(labelsize=19)         
plt.plot(tlist, zs, marker = '.', color = 'black', label = 'Exact',linewidth=4)
plt.plot(t1[step1],mz1[step1], marker = '*', markersize = 10, color = 'blue',label = 'Ansatz II Rank=1',linewidth=4)
plt.plot(t8[step1],mz8[step1], marker = 'o', markersize = 10, color = 'red',label = 'Ansatz II Rank=8',linewidth=4)

plt.legend(prop = font1)

plt.title('Average magnetization of TFIM ',fontsize=24)
plt.xlabel('t', fontsize = 24)
plt.ylabel('Average magnetization', fontsize = 24)
figname = 'mzTFIM.png'
plt.savefig(figname)

plt.figure(figsize = (12.0,9.0))
plt.tick_params(labelsize=22)         
plt.plot(tlist, xs, marker = '.', color = 'black', label = 'Exact',linewidth=4)
plt.plot(t1[step1],mx1[step1], marker = '*', markersize = 12, color = 'blue',label = 'Ansatz II Rank=1',linewidth=4)
plt.plot(t8[step1],mx8[step1], marker = '.', markersize = 10, color = 'red',label = 'Ansatz II Rank=8',linewidth=4)
plt.legend(prop = font1)

# plt.plot(gamma*tlist, rho2, marker='.')
plt.title('Average <$\sigma^x$> of TFIM ',fontsize=24)
plt.xlabel('t', fontsize = 24)
plt.ylabel('Average <$\sigma^x$>', fontsize = 24)
figname = 'mxTFIM.png'
plt.savefig(figname)
