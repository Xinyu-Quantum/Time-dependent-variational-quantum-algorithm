# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 21:22:14 2022

@author: xinyu
"""

import numpy as np
from matplotlib import pyplot as plt
import pickle


dataname = 'TLSqutip.pkl'
f = open(dataname,'rb')
h,gamma,spin_num,medatas,tlist,excitation_num,rho2q = pickle.load(f)
f.close()

dataname = 'TLSexactver2.pkl'
f = open(dataname,'rb')
nev2, rho2v2, tv2, paras, mats, vecs = pickle.load(f)
f.close()

dataname = 'TLSexact.pkl'
f = open(dataname,'rb')
ne, rho2, rho01, tv1, paras, mats, vecs = pickle.load(f)
f.close()

datanamever1 = 'TLSny'
ne_ny1 = [[] for i in range(10)]
rho2_ny1 = [[] for i in range(10)] 
for i in range (10):
    dataname = datanamever1+str(i+1)+'.pkl'
    f = open(dataname,'rb')
    ne_ny1[i], rho2_ny1[i], t, paras, mats, vecs = pickle.load(f)
    f.close()
ne_ny1 = np.array(ne_ny1)
ne_ny_1 = np.mean(ne_ny1,axis = 0)
ne_ny_1_std = np.std(ne_ny1,axis = 0)
rho2_ny1 = np.array(rho2_ny1)
rho2_ny_1 = np.mean(rho2_ny1,axis = 0)
rho2_ny_1_std = np.std(rho2_ny1,axis = 0)

datanamever2 = 'TLSnyver2'
ne_ny2 = [[] for i in range(10)]
rho2_ny2 = [[] for i in range(10)] 
for i in range (10):
    dataname = datanamever2+str(i+1)+'.pkl'
    f = open(dataname,'rb')
    ne_ny2[i], rho2_ny2[i], t, paras, mats, vecs = pickle.load(f)
    f.close()
ne_ny2 = np.array(ne_ny2)
ne_ny_2 = np.mean(ne_ny2,axis = 0)
ne_ny_2_std = np.std(ne_ny2,axis = 0)
rho2_ny2 = np.array(rho2_ny2)
rho2_ny_2 = np.mean(rho2_ny2,axis = 0)
rho2_ny_2_std = np.std(rho2_ny2,axis = 0)

font1 = {
    'family':'Times New Roman',
    'weight': 'normal',
    'size': 24,
    }

step1 = np.linspace(0,490,50, dtype = int)
step2 = np.linspace(0,192,25, dtype = int)
step3 = np.linspace(0,240,25,  dtype = int)

###
plt.figure(figsize = (12.0,9.0))
plt.tick_params(labelsize=22)         
plt.plot(tlist, excitation_num, marker='.', color = 'black', label = 'Exact')
plt.plot(tv1[step1],ne[step1],marker = '*', markersize = 10, color = 'red', label = 'Ansatz I SV')
plt.plot(tv2[step1],nev2[step1], marker = 'o', linestyle = '', markersize = 10, color = 'blue', label = 'Ansatz II SV')
plt.legend(prop = font1, loc = 'center right')

# plt.plot(gamma*tlist, rho2, marker='.')
plt.xlabel('t', fontsize = 24)
plt.ylabel('Excitation population', fontsize = 24)
plt.title('Excitation population of TLS', fontsize = 24)
figname = 'TLSexc.png'
plt.savefig(figname)

###
plt.figure(figsize = (12.0,9.0))
plt.tick_params(labelsize=22)         
plt.plot(tlist, rho2q, marker='.', color = 'black', label = 'Exact')
plt.plot(tv1[step1],rho2[step1],marker = '*', markersize = 10, color = 'red', label = 'Ansatz I SV')
plt.plot(tv2[step1],rho2v2[step1], marker = 'o', linestyle = '', markersize = 10, color = 'blue', label = 'Ansatz II SV')
plt.legend(prop = font1, loc = 'center right')


plt.xlabel('t', fontsize = 24)
plt.ylabel('Tr(ρ$^2$)', fontsize = 24)
plt.title('Tr[ρ$^2$] of TLS', fontsize = 24)
figname = 'TLSrho2.png'
plt.savefig(figname)

###
plt.figure(figsize = (12.0,9.0))
plt.tick_params(labelsize=22)         
plt.plot(tlist, excitation_num, marker='.', color = 'black', label = 'Exact')
plt.errorbar(t[step1], ne_ny_1[step1], yerr = ne_ny_1_std[step1], uplims=True, lolims=True, marker = '*', color = 'red', ecolor = 'red', label = 'Ansatz I NM')
plt.errorbar(t[step1], ne_ny_2[step1], yerr = ne_ny_2_std[step1], uplims=True, lolims=True, marker = 'o', color = 'blue', ecolor = 'blue', label = 'Ansatz II NM')
plt.legend(prop = font1, loc = 'center right')


plt.xlabel('t', fontsize = 24)
plt.ylabel('Excitation population', fontsize = 24)
plt.title('Excitation population of TLS', fontsize = 24)
figname = 'TLSexcny.png'
plt.savefig(figname)


###
plt.figure(figsize = (12.0,9.0))
plt.tick_params(labelsize=22)         
plt.plot(tlist, rho2q, marker='.', color = 'black', label = 'Exact')
plt.errorbar(t[step1], rho2_ny_1[step1], yerr = rho2_ny_1_std[step1], uplims=True, lolims=True, marker = '*', color = 'red', ecolor = 'red', label = 'Ansatz I NM')
plt.errorbar(t[step1], rho2_ny_2[step1], yerr = rho2_ny_2_std[step1], uplims=True, lolims=True, marker = 'o', color = 'blue', ecolor = 'blue', label = 'Ansatz II NM')
plt.legend(prop = font1, loc = 'center right')


plt.xlabel('t', fontsize = 24)
plt.ylabel('Tr(ρ$^2$)', fontsize = 24)
plt.title('Tr[ρ$^2$] of TLS', fontsize = 24)
figname = 'TLSrho2ny.png'
plt.savefig(figname)
