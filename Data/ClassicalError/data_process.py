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
h,gamma,spin_num,Rs,Ents,tlist,Mats,zs,xs = pickle.load(f)
f.close()

data1 = loadmat('t0.02TA_N2L2R4Jz1h1.mat')

a = np.linspace(0,249,250,dtype = int)

t1 = data1['tlist']
t1 = t1[0][a]
mz1 = data1['mz']
mz1 = mz1[0][a]
mx1 = data1['mx']
mx1 = mx1[0][a]

a = np.linspace(0,499,500,dtype = int)
data2 = loadmat('t0.01TA_N2L2R4Jz1h1.mat')
t2 = data2['tlist']
t2 = t2[0][a]
mz2 = data2['mz']
mz2 = mz2[0][a]
mx2 = data2['mx']
mx2 = mx2[0][a]

a = np.linspace(0,4999,5000,dtype = int)
data3 = loadmat('t0.001TA_N2L2R4Jz1h1.mat')
t3 = data3['tlist']
t3 = t3[0][a]
mz3 = data3['mz']
mz3 = mz3[0][a]
mx3 = data3['mx']
mx3 = mx3[0][a]

step1 = np.linspace(0,490,50, dtype = int)

font1 = {
    'family':'Times New Roman',
    'weight': 'normal',
    'size': 24,
    }

plt.figure(figsize = (12.0,9.0))
plt.tick_params(labelsize=22)         
plt.plot(tlist, zs, marker = '.', color = 'black', label = 'Exact',linewidth=4)
plt.plot(t1,mz1, marker = '.', color = 'blue',label = 'Ansatz II dt=0.02',linewidth=4)
plt.plot(t2,mz2, marker = '.', color = 'red',label = 'Ansatz II dt=0.01',linewidth=4)
plt.plot(t3,mz3, marker = '.', color = 'green',label = 'Ansatz II dt=0.001',linewidth=4)

plt.legend(prop = font1)

plt.title('Average magnetization of TFIM ',fontsize=24)
plt.xlabel('t', fontsize = 24)
plt.ylabel('Average magnetization', fontsize = 24)
figname = 'mzdtTFIM.png'
plt.savefig(figname)

plt.figure(figsize = (12.0,9.0))
plt.tick_params(labelsize=22)         
plt.plot(tlist, xs, marker = '.', color = 'black', label = 'Exact',linewidth=4)
plt.plot(t1,mx1, marker = '.', color = 'blue',label = 'Ansatz II dt=0.02',linewidth=4)
plt.plot(t2,mx2, marker = '.', color = 'red',label = 'Ansatz II dt=0.01',linewidth=4)
plt.plot(t3,mx3, marker = '.', color = 'green',label = 'Ansatz II dt=0.001',linewidth=4)

plt.legend(prop = font1)

# plt.plot(gamma*tlist, rho2, marker='.')
plt.title('Average <$\sigma^x$> of TFIM ',fontsize=24)
plt.xlabel('t', fontsize = 24)
plt.ylabel('Average <$\sigma^x$>', fontsize = 24)
figname = 'mxdtTFIM.png'
plt.savefig(figname)
