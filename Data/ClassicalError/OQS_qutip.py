# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 13:33:04 2022

@author: xinyu
"""

from qutip import*
import numpy as np
from matplotlib import pyplot as plt
import pickle

def generate_Hamitonian_XYZ(spin_num, Jx, Jy, Jz ,h):
    Hxx = [[] for i in range(spin_num)]
    Hyy = [[] for i in range(spin_num)]
    Hzz = [[] for i in range(spin_num)]    
    for i in range(spin_num):
        if(i == spin_num-1):
            Hxx[-1] = tensor(sigmax(),qeye(2**(spin_num-2)))
            Hxx[-1] = tensor(Hxx[-1],sigmax())
            
            Hyy[-1] = tensor(sigmay(),qeye(2**(spin_num-2)))
            Hyy[-1] = tensor(Hyy[-1],sigmay())
            
            Hzz[-1] = tensor(sigmaz(),qeye(2**(spin_num-2)))
            Hzz[-1] = tensor(Hzz[-1],sigmaz())
        else:
            Hxx[i] = tensor(qeye(2**i),sigmax())
            Hxx[i] = tensor(Hxx[i],sigmax())
            Hxx[i] = tensor(Hxx[i],qeye(2**(spin_num-i-2)))
        
            Hyy[i] = tensor(qeye(2**(i)),sigmay())
            Hyy[i] = tensor(Hyy[i],sigmay())
            Hyy[i] = tensor(Hyy[i],qeye(2**(spin_num-i-2)))
        
            Hzz[i] = tensor(qeye(2**(i)),sigmaz())
            Hzz[i] = tensor(Hzz[i],sigmaz())
            Hzz[i] = tensor(Hzz[i],qeye(2**(spin_num-i-2)))
        #print(Hzz[i])
        # if(i%spin_num!=0 and (i+1)%spin_num!=0):
        #     Hxx[i] = qeye(2)
        #     Hyy[i] = qeye(2)
        #     Hzz[i] = qeye(2)
        # else:
        #     Hxx[i] = sigmax()
        #     Hyy[i] = sigmay()
        #     Hzz[i] = sigmaz()
    # for j in range(spin_num-1):
    #     for i in range(spin_num):
    #         if(i%spin_num!=j and (i+1)%spin_num!=j):
    #             Hxx[i] = tensor(Hxx[i],qeye(2))
    #             Hyy[i] = tensor(Hyy[i],qeye(2))
    #             Hzz[i] = tensor(Hzz[i],qeye(2))
    #         else:
    #             Hxx[i] = tensor(Hxx[i],sigmax())
    #             Hyy[i] = tensor(Hyy[i],sigmay())
    #             Hzz[i] = tensor(Hzz[i],sigmaz())
    
              
    Hx = [[] for i in range(spin_num)]
    for i in range(spin_num):
        Hx[i] = tensor(qeye(2**i),sigmax())
        Hx[i] = tensor(Hx[i],qeye(2**(spin_num-i-1)))
        #print(Hx[i])
    # Hx[0] = sigmax()
    # for i in range(spin_num-1):
    #     Hx[i+1] = qeye(2)
    # for i in range(spin_num-1):
    #     for j in range(spin_num):
    #         if(j == i+1):
    #             Hx[j] = tensor(Hx[j],sigmax())
    #         else:
    #             Hx[j] = tensor(Hx[j],qeye(2))
                
    H = 0 
        
    if(spin_num == 2):
        H += Jx*Hxx[0].full()+Jy*Hyy[0].full()+Jz*Hzz[0].full()+h*Hx[0].full()+h*Hx[1].full()
    elif(spin_num>2):
        for i in range(spin_num):
            H += Jx*Hxx[i].full()+Jy*Hyy[i].full()+Jz*Hzz[i].full()+h*Hx[i].full()         
            #H = H-Jz*Hzz[-1].full()
    return Qobj(H)

def generate_cops(spin_num, gamma):
    cops = [[] for i in range(spin_num)]
    for i in range(spin_num):
        cops[i] = tensor(qeye(2**i),sigmam())
        cops[i] = tensor(cops[i],qeye(2**(spin_num-i-1)))
    # for i in range(spin_num-1):
    #     cops[i+1] = qeye(2)
    # cops[0] = sigmam()
    # for i in range(spin_num-1):
    #     for j in range(spin_num):
    #         if(j == i+1):
    #             cops[j] = tensor(cops[j],sigmam())
    #         else:
    #             cops[j] = tensor(cops[j],qeye(2))
    c_ops = []
    for i in range(spin_num):
        c_ops.append(Qobj(np.sqrt(gamma)*cops[i].full())) 
    return c_ops
  
def generate_eops_sigmaz(spin_num):
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
    e_ops.append(Qobj(e_op.full()/spin_num))
    return e_ops
    
def generate_eops_sigmax(spin_num):
    eops = [[] for i in range(spin_num)]
    for i in range(spin_num-1):
        eops[i+1] = qeye(2)
    eops[0] = sigmax()
    for i in range(spin_num-1):
        for j in range(spin_num):
            if(j == i+1):
                eops[j] = tensor(eops[j],sigmax())
            else:
                eops[j] = tensor(eops[j],qeye(2))
    e_ops = []
    e_op = 0
    for i in range(spin_num):
        e_op += eops[i]
    e_ops.append(Qobj(e_op.full()/spin_num))
    return e_ops

def generate_eops_sigmay(spin_num):
    eops = [[] for i in range(spin_num)]
    for i in range(spin_num-1):
        eops[i+1] = qeye(2)
    eops[0] = sigmay()
    for i in range(spin_num-1):
        for j in range(spin_num):
            if(j == i+1):
                eops[j] = tensor(eops[j],sigmay())
            else:
                eops[j] = tensor(eops[j],qeye(2))
    e_ops = []
    e_op = 0
    for i in range(spin_num):
        e_op += eops[i]
    e_ops.append(Qobj(e_op.full()/spin_num))
    return e_ops

def generate_psi0(spin_num):
    psi0 = basis(2,0)
    for i in range(spin_num-1):
        psi0 = tensor(psi0,basis(2,0))
    return Qobj(psi0.full())

def generate_psi1(spin_num):
    psi1 = basis(2,1)
    for i in range(spin_num-1):
        psi1 = tensor(psi1,basis(2,1))
    return Qobj(psi1.full())

def Ent_cal(rho_diagonal):
    rho_diagonal = rho_diagonal[rho_diagonal>0]
    return np.sum(-rho_diagonal*np.log(rho_diagonal))

def Rank_rho(rho_list, spin_num):
    R = np.zeros(len(rho_list))
    Ent = np.zeros(len(rho_list))
    Mats = np.zeros([len(rho_list),2**spin_num,2**spin_num],dtype = complex)
    for i in range(len(rho_list)):
        Mats[i] = rho_list[i].full(order='C')
        w,v = np.linalg.eig(rho_list[i].full(order='C'))
        w = np.real(w)
        print(w)
        Ent[i] = Ent_cal(w)/spin_num
        R[i] = np.count_nonzero(w>0.001)
    return Mats, R, Ent

Jx = 0
Jy = 0
Jz = 1.
h = 1
gamma = 0.1
spin_num = 2


tlist = np.linspace(0., 5., 500)
H = generate_Hamitonian_XYZ(spin_num, Jx, Jy, Jz, h)
c_ops = generate_cops(spin_num,gamma)
psi1 = generate_psi1(spin_num)
psi0 = generate_psi0(spin_num)
medatas = mesolve(H, psi0, tlist, c_ops, [])  
Mats, Rs, Ents = Rank_rho(medatas.states, spin_num)
e_ops_z = generate_eops_sigmaz(spin_num)
e_ops_x = generate_eops_sigmax(spin_num)
e_ops_y = generate_eops_sigmay(spin_num)
zs = expect(e_ops_z, medatas.states)[0]
xs = expect(e_ops_x, medatas.states)[0]
ys = expect(e_ops_y, medatas.states)[0]

plt.figure(figsize = (12.0,9.0))
plt.tick_params(labelsize=20)         
plt.plot(tlist, zs, marker='.')
plt.title('Average magnetization of TFIM ',fontsize=22)
plt.xlabel('t', fontsize = 22)
plt.ylabel('Average magnetization', fontsize = 22)
figname = 'TFIMqutip.jpg'
plt.savefig(figname)

plt.figure(figsize = (12.0,9.0))
plt.tick_params(labelsize=20)         
plt.plot(tlist, xs, marker='.')
plt.title('Average <$\sigma^x$> of TFIM ',fontsize=22)
plt.xlabel('t', fontsize = 22)
plt.ylabel('Average <$\sigma^x$>', fontsize = 22)
figname = 'mxTFIMqutip.jpg'
plt.savefig(figname)

plt.figure(figsize = (12.0,9.0))
plt.tick_params(labelsize=20)         
plt.plot(tlist, xs, marker='.')
plt.title('Average <$\sigma^y$> of TFIM ',fontsize=22)
plt.xlabel('t', fontsize = 22)
plt.ylabel('Average <$\sigma^y$>', fontsize = 22)
figname = 'myTFIMqutip.jpg'
plt.savefig(figname)

plt.figure(figsize = (12.0,9.0))
plt.tick_params(labelsize=20)         
plt.plot(gamma*tlist, Rs, marker='.')
plt.title('h'+str(h)+'gamma'+str(gamma)+'spin'+str(spin_num),fontsize=22)
plt.xlabel('γt', fontsize = 22)
plt.ylabel('Rank', fontsize = 22)
figname = 'h'+str(h)+'gamma'+str(gamma)+'spin'+str(spin_num)+'Rank.jpg'
plt.savefig(figname)
    
plt.figure(figsize = (12.0,9.0))
plt.tick_params(labelsize=20)         
plt.plot(gamma*tlist, Ents, marker='.')
plt.title('h'+str(h)+'gamma'+str(gamma)+'spin'+str(spin_num),fontsize=22)
plt.xlabel('γt', fontsize = 22)
plt.ylabel('Entropy', fontsize = 22)
figname = 'h'+str(h)+'gamma'+str(gamma)+'spin'+str(spin_num)+'Entropy.jpg'
plt.savefig(figname)
    
dataname = 'TFIMqutip.pkl'
f = open(dataname,'wb')
pickle.dump((h,gamma,spin_num,Rs,Ents,tlist,Mats,zs,xs),f)
f.close() 