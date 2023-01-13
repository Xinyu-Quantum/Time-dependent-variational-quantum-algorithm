# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 21:22:14 2022

@author: xinyu
"""

import numpy as np
from matplotlib import pyplot as plt
import pickle


dataname = 'TFIMqutip.pkl'
f = open(dataname,'rb')
h,gamma,spin_num,tlist,zs = pickle.load(f)
f.close()

dataname = 'TFIMexact.pkl'
f = open(dataname,'rb')
mz, t, paras, mats, vecs = pickle.load(f)
f.close()

dataname = 'TFIMexactver2.pkl'
f = open(dataname,'rb')
mzv2, tv2, parasv2, matsv2, vecsv2 = pickle.load(f)
f.close()

# dataname = '0.02RegTFIMexact.pkl'
# f = open(dataname,'rb')
# mz3, t, paras, mats, vecs = pickle.load(f)
# f.close()

# dataname = '0.02RegTFIMexactver2.pkl'
# f = open(dataname,'rb')
# mz4, tv2, parasv2, matsv2, vecsv2 = pickle.load(f)
# f.close()

font1 = {
    'family':'Times New Roman',
    'weight': 'normal',
    'size': 24,
    }

step1 = np.linspace(0,490,50, dtype = int)

datanamever1 = 'TFIMny'
mz_ny1 = [[] for i in range(10)]
for i in range (10):
    dataname = datanamever1+str(i+1)+'.pkl'
    f = open(dataname,'rb')
    mz_ny1[i], t, paras, mats, vecs = pickle.load(f)
    f.close()
mz_ny1 = np.array(mz_ny1)
mz_ny_1 = np.mean(mz_ny1,axis = 0)
mz_ny_1_std = np.std(mz_ny1,axis = 0)

datanamever2 = 'TFIMver2.pkl'
mz_ny2 = [[] for i in range(7)]
for i in range (7):
    dataname = 'Ny'+str(i+1)+datanamever2
    f = open(dataname,'rb')
    mz_ny2[i], t, paras, mats, vecs = pickle.load(f)
    f.close()
mz_ny2 = np.array(mz_ny2)
mz_ny_2 = np.mean(mz_ny2,axis = 0)
mz_ny_2_std = np.std(mz_ny2,axis = 0)


plt.figure(figsize = (12.0,9.0))
plt.tick_params(labelsize=22)         
plt.plot(tlist, zs, marker='.', color = 'black', label = 'Exact')
plt.plot(t,mz,marker = '.', color = 'red', label = 'Ansatz I SV')
plt.plot(tv2,mzv2, marker = '.', color = 'blue', label = 'Ansatz II SV')
plt.legend(prop = font1, loc = 'upper right')

plt.title('Average magnetization of TFIM ',fontsize=24)
plt.xlabel('t', fontsize = 24)
plt.ylabel('Average magnetization', fontsize = 24)
figname = 'TFIM.png'
plt.savefig(figname)

# plt.figure(figsize = (12.0,9.0))
# plt.tick_params(labelsize=22)         
# plt.plot(tlist, zs, marker='.', color = 'black', label = 'Exact')
# plt.plot(t,mz3,marker = '.', color = 'red', label = 'Ansatz I SV Reg=0.02')
# plt.plot(tv2,mz4, marker = '.', color = 'blue', label = 'Ansatz II SV Reg=0.02')
# plt.legend(prop = font1, loc = 'upper right')

# plt.title('Average magnetization of TFIM ',fontsize=24)
# plt.xlabel('t', fontsize = 24)
# plt.ylabel('Average magnetization', fontsize = 24)
# figname = 'RegTFIM.png'
# plt.savefig(figname)

plt.figure(figsize = (12.0,9.0))
plt.tick_params(labelsize=22)         
plt.plot(tlist, zs, marker='.', color = 'black', label = 'Exact')
plt.errorbar(t[step1],mz_ny_1[step1], yerr= mz_ny_1_std[step1], uplims=True, lolims=True, marker = '*', color = 'red', ecolor = 'red', label = 'Ansatz I NM')
plt.errorbar(t[step1],mz_ny_2[step1], yerr= mz_ny_2_std[step1], uplims=True, lolims=True, marker = 'o', color = 'blue', ecolor = 'blue', label = 'Ansatz II NM')
#plt.plot(tv2,mzv2, marker = 'o', color = 'blue', label = 'Ansatz II SV')
plt.legend(prop = font1, loc = 'upper right')

plt.title('Average magnetization of TFIM ',fontsize=24)
plt.xlabel('t', fontsize = 24)
plt.ylabel('Average magnetization', fontsize = 24)
figname = 'TFIMny.png'
plt.savefig(figname)

