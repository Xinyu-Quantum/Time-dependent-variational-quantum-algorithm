# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 13:33:04 2022

@author: xinyu
"""

from qutip import*
import numpy as np
from matplotlib import pyplot as plt
import pickle

def generate_Hamitonian(spin_num, J ,h):
    H0 = [[] for i in range(spin_num-1)]
    if (spin_num>2):
        H0[0] = tensor(sigmaz(),sigmaz())
        H0[1] = tensor(qeye(2),sigmaz()) 
        for i in range(spin_num-2):
            H0[i+2] = tensor(qeye(2),qeye(2))
            for i in range(spin_num-2):
                for j in range(spin_num):
                    if(j == i+1):
                        H0[j] = tensor(H0[j],sigmaz())
                    elif(j == i+2):
                        H0[j] = tensor(H0[j],sigmaz())
                    else:
                        H0[j] = tensor(H0[j],qeye(2))
    
    if(spin_num == 2):
        H0[0] = tensor(sigmaz(),sigmaz())
    
    H1 = [[] for i in range(spin_num)]
    H1[0] = sigmax()
    for i in range(spin_num-1):
        H1[i+1] = qeye(2)
    for i in range(spin_num-1):
        for j in range(spin_num):
            if(j == i+1):
                H1[j] = tensor(H1[j],sigmax())
            else:
                H1[j] = tensor(H1[j],qeye(2))
                
    H = 0
    for i in range(spin_num-1):
        H += -J*H0[i]-h*H1[i]
    H += -h*H1[-1]
    return H
    
def generate_cops(spin_num, gamma):
    cops = [[] for i in range(spin_num)]
    for i in range(spin_num-1):
        cops[i+1] = qeye(2)
    cops[0] = sigmam()
    for i in range(spin_num-1):
        for j in range(spin_num):
            if(j == i+1):
                cops[j] = tensor(cops[j],sigmam())
            else:
                cops[j] = tensor(cops[j],qeye(2))
    c_ops = []
    for i in range(spin_num):
        c_ops.append(np.sqrt(gamma)*cops[i])
    return c_ops
  
def generate_eops(spin_num):
    eops = [[] for i in range(spin_num)]
    for i in range(spin_num-1):
        eops[i+1] = qeye(2)
    eops[0] = sigmaz()
    for i in range(spin_num-1):
        for j in range(spin_num):
            if(j == i+1):
                eops[j] = tensor(eops[j],sigmaz())
            else:
                eops[j] = tensor(eops[j],qeye(2))
    e_ops = []
    e_op = 0
    for i in range(spin_num):
        e_op += eops[i]
    e_ops.append(e_op/spin_num)
    return e_ops
    
def generate_psi0(spin_num):
    psi0 = basis(2,0)
    for i in range(spin_num-1):
        psi0 = tensor(psi0,basis(2,0))
    return psi0

def generate_psi1(spin_num):
    psi1 = basis(2,1)
    for i in range(spin_num-1):
        psi1 = tensor(psi1,basis(2,1))
    return psi1

def Ent_cal(rho_diagonal):
    rho_diagonal = rho_diagonal[rho_diagonal!=0]
    return np.sum(-rho_diagonal*np.log(rho_diagonal))

def Rank_rho(rho_list, spin_num):
    R = np.zeros(len(rho_list))
    Ent = np.zeros(len(rho_list))
    for i in range(len(rho_list)):
        w,v = np.linalg.eig(rho_list[i].full(order='C'))
        w = np.real(w)
        Ent[i] = Ent_cal(w)/spin_num
        R[i] = np.count_nonzero(w>1./2**spin_num)
    return R, Ent


J = 1.
h = 1.
gamma = 0.1
spin_num = 2
e_ops = generate_eops(spin_num)

tlist = np.linspace(0., 5., 500)
H = generate_Hamitonian(spin_num, J, h)
c_ops = generate_cops(spin_num,gamma)
psi0 = generate_psi0(spin_num)
medatas = mesolve(H, psi0, tlist, c_ops, e_ops)  

zs = medatas.expect[0]
    
plt.figure(figsize = (12.0,9.0))
plt.tick_params(labelsize=20)         
plt.plot(tlist, zs, marker='.')
plt.title('Average magnetization of TFIM ',fontsize=22)
plt.xlabel('t', fontsize = 22)
plt.ylabel('Average magnetization', fontsize = 22)
figname = 'TFIMqutip.jpg'
plt.savefig(figname)
    

    
dataname = 'TFIMqutip.pkl'
f = open(dataname,'wb')
pickle.dump((h,gamma,spin_num,tlist,zs),f)
f.close() 