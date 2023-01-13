# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 16:36:16 2022

@author: xinyu
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 14:38:54 2022

@author: xinyu
"""

import numpy as np
from math import pi
import pickle
import time
from matplotlib import pyplot as plt
from itertools import product

from qiskit import*
from qiskit.compiler import transpile, assemble
from qiskit.tools.jupyter import *
from qiskit.visualization import *
from qiskit.circuit import Parameter
from qiskit.providers.aer.noise import NoiseModel
from qiskit.providers.aer.noise.errors.standard_errors import depolarizing_error, coherent_unitary_error, amplitude_damping_error, phase_damping_error
from qiskit import IBMQ
from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister, execute
from qiskit.transpiler.passes import ALAPSchedule, DynamicalDecoupling, RemoveBarriers, RemoveFinalMeasurements

def obtain_counts(counts, qstr):
    value = counts.get(qstr)
    if(value == None):
        value = 0
    return value

IBMQ.save_account("7f647333d63485cd69118266daf5fabd1eb01ccc82ddca48ab094afbe21669a69db036d6dd5bde3e4f42ef6d2aa343280e87ca9389bdaed6c587ac1f9bfffd1b");

#import noise model
provider = IBMQ.load_account()


# define backend
qasm_simulator = Aer.get_backend('qasm_simulator')
statevector_simulator = Aer.get_backend('statevector_simulator')
noisy_simulator = provider.get_backend('ibmq_manila')

manila_noise = NoiseModel.from_backend(noisy_simulator)
basis_gates = manila_noise.basis_gates

qubit_num = 2

diagpara_num = 1
vecpara_num = 2

Parameters = [Parameter('γ'),Parameter('θ0'),Parameter('θ1'),Parameter('θ2'),Parameter('θ3')]
para = np.array([0.,0.,0.,0.,0.])
ele_num = [2,4,4,4,4]

circuit_elements =  [[]for i in range(len(Parameters))]
circuit_elements[0] = QuantumCircuit(2*qubit_num,name = 'ry')
circuit_elements[0].ry(Parameters[0],[0,1])
circuit_elements[0].cx(0,2)
circuit_elements[0].cx(1,3)
circuit_elements[1] = QuantumCircuit(2*qubit_num,name = 'rx0')
circuit_elements[1].rx(Parameters[1],[0,1])
circuit_elements[1].rx(-1*Parameters[1],[2,3])
circuit_elements[2] = QuantumCircuit(2*qubit_num,name = 'ry0')
circuit_elements[2].ry(Parameters[2],[0,1])
circuit_elements[2].ry(Parameters[2],[2,3])
circuit_elements[2].cz(0,1)
circuit_elements[2].cz(2,3)
circuit_elements[3] = QuantumCircuit(2*qubit_num,name = 'rx1')
circuit_elements[3].rx(Parameters[3],[0,1])
circuit_elements[3].rx(-1*Parameters[3],[2,3])
circuit_elements[4] = QuantumCircuit(2*qubit_num,name = 'ry1')
circuit_elements[4].ry(Parameters[4],[0,1])
circuit_elements[4].ry(Parameters[4],[2,3])
circuit_elements[4].cz(0,1)
circuit_elements[4].cz(2,3)

# circuit_elements[3] = QuantumCircuit(2*qubit_num,name = 'rx1')
# circuit_elements[3].rzz(Parameters[3],0,1)
# circuit_elements[3].rzz(-1*Parameters[3],2,3)

diff_elements = [[]for i in range(len(Parameters))]
diff_elements[0] = [[]for i in range(ele_num[0])]
for i in range(ele_num[0]):
    diff_elements[0][i] = QuantumCircuit(2*qubit_num+1,name = 'cy0'+str(i+1))
    diff_elements[0][i].cy(0,i+1)

diff_elements[1] =  [[]for i in range(ele_num[1])]
for i in range(ele_num[1]):
    diff_elements[1][i] = QuantumCircuit(2*qubit_num+1,name = 'cx0'+str(i+1))
    diff_elements[1][i].cx(0,i+1)
    
diff_elements[2] =  [[]for i in range(ele_num[2])]
for i in range(ele_num[2]):
    diff_elements[2][i] = QuantumCircuit(2*qubit_num+1,name = 'cy0'+str(i+1))
    diff_elements[2][i].cy(0,i+1)
    
diff_elements[3] =  [[]for i in range(ele_num[3])]
for i in range(ele_num[3]):
    diff_elements[3][i] = QuantumCircuit(2*qubit_num+1,name = 'cx0'+str(i+1))
    diff_elements[3][i].cx(0,i+1)

diff_elements[4] =  [[]for i in range(ele_num[4])]
for i in range(ele_num[4]):
    diff_elements[4][i] = QuantumCircuit(2*qubit_num+1,name = 'cy0'+str(i+1))
    diff_elements[4][i].cy(0,i+1)
    
# diff_elements[3] =  [[]for i in range(ele_num[3])]
# for i in range(ele_num[3]):
#     diff_elements[3][i] = QuantumCircuit(2*qubit_num+1,name = 'czz0')
#     diff_elements[3][i].cz(0,2*i+1)
#     diff_elements[3][i].cz(0,2*i+2)
    
    
diff_coeff = [[]for i in range(len(Parameters))]
diff_coeff[0] = [[]for i in range(ele_num[0])]
for i in range(ele_num[0]):
    diff_coeff[0][i] = 1
    
diff_coeff[1] = [[]for i in range(ele_num[1])]
for i in range(2):
    diff_coeff[1][i] = 1
    diff_coeff[1][i+2] = -1
        
diff_coeff[2] = [[]for i in range(ele_num[2])]
for i in range(ele_num[2]):
    diff_coeff[2][i] = 1
    
diff_coeff[3] = [[]for i in range(ele_num[3])]
for i in range(2):
    diff_coeff[3][i] = 1
    diff_coeff[3][i+2] = -1
    
diff_coeff[4] = [[]for i in range(ele_num[4])]
for i in range(ele_num[4]):
    diff_coeff[4][i] = 1
    
# diff_coeff[3] = [[]for i in range(ele_num[3])]
# diff_coeff[3][0] = 1
# diff_coeff[3][1] = -1

# in the case of our circuit, the product of two coefficient is always -0.25
diff_coefficient_prod = 0.25

def generate_circuit_matrix(index1, index2, ele1, ele2):
    circ = []
    if(index1<index2):
        ancilla = QuantumRegister(1)
        phys = QuantumRegister(2*qubit_num)
        circ = QuantumCircuit(ancilla,phys)
        
        circ.h(0)
        circ.x(0)
        if(index1>0):
            for i in range(index1):
                circ.append(circuit_elements[i].to_instruction(),[phys[i] for i in range(2*qubit_num)])
        circ.barrier()
        circ.append(diff_elements[index1][ele1].to_instruction(),list(range(2*qubit_num+1)))
        circ.barrier()
        circ.x(0)
        for i in range(index1,index2):
            circ.append(circuit_elements[i].to_instruction(),[phys[i] for i in range(2*qubit_num)])
        circ.barrier()
        circ.append(diff_elements[index2][ele2].to_instruction(),list(range(2*qubit_num+1)))
        circ.barrier()
        circ.h(0)
        
    if(index1 == index2):
        ancilla = QuantumRegister(1)
        phys = QuantumRegister(2*qubit_num)
        circ = QuantumCircuit(ancilla,phys)
        
        circ.h(0)
        circ.x(0)
        if(index1>0):
            for i in range(index1):
                circ.append(circuit_elements[i].to_instruction(),[phys[i] for i in range(2*qubit_num)])
        circ.barrier()
        circ.append(diff_elements[index1][ele1].to_instruction(),list(range(2*qubit_num+1)))
        circ.barrier()
        circ.x(0)
        circ.barrier()
        circ.append(diff_elements[index2][ele2].to_instruction(),list(range(2*qubit_num+1)))
        circ.barrier()
        circ.h(0)
           
    return circ

circuit_matrix = [[] for i in range(len(Parameters))]
coeff_matrix = [[] for i in range(len(Parameters))]
for i in range(len(Parameters)):
    circuit_matrix[i] = [[]for j in range(i+1)]
    coeff_matrix[i] = [[]for j in range(i+1)]
        
for i in range(len(Parameters)):
    for j in range(i+1):
        circuit_matrix[i][j] = []
        coeff_matrix[i][j] = []
        for p,q in product(range(ele_num[i]),range(ele_num[j])):
            circuit_matrix[i][j].append(generate_circuit_matrix(j,i,q,p))
            coeff_matrix[i][j].append(diff_coeff[j][q]*diff_coeff[i][p])

def evaluate_matrix(circuit_matrix, para, coeff_matrix, shot_num):
    mat = np.ones([len(para),len(para)],dtype = 'float')
    for i in range(len(Parameters)):
        for j in range(i+1):    
            value = 0
            for k, circ in enumerate(circuit_matrix[i][j]):
                circuit = circ.bind_parameters({Parameters[m]: para[m] for m in range(i)})
                circuit = transpile(circuit, basis_gates = basis_gates,optimization_level=3)
                job = execute(circuit,backend = statevector_simulator)     
                result = job.result().get_statevector()
                p = 0
                for m in range(2**(2*qubit_num)):
                    p += abs(result[2*m])**2
                value += (2*p-1)* coeff_matrix[i][j][k]
            # if(abs(value)<10**(-4)):
            #     value = 10**(-4)
            mat[i][j] = value
            mat[j][i] = value
    return mat*0.25

mat = evaluate_matrix(circuit_matrix,para,coeff_matrix,2**16)

J = 1
h = 1
gamma = 0.1

vec_elements = []
vec_coefficient = []
vec_phase = []

#0
vecirc = QuantumCircuit(2*qubit_num+1)
vecirc.cz(0,1)
vecirc.cz(0,2)
vec_elements.append(vecirc)
vec_coefficient.append(-0.5*J)
vec_phase.append(0)

#1
vecirc = QuantumCircuit(2*qubit_num+1)
vecirc.cz(0,3)
vecirc.cz(0,4)
vec_elements.append(vecirc)
vec_coefficient.append(0.5*J)
vec_phase.append(0)

#2
vecirc = QuantumCircuit(2*qubit_num+1)
vecirc.cx(0,1)
vec_elements.append(vecirc)
vec_coefficient.append(-0.5*h)
vec_phase.append(0)

#3
vecirc = QuantumCircuit(2*qubit_num+1)
vecirc.cx(0,2)
vec_elements.append(vecirc)
vec_coefficient.append(-0.5*h)
vec_phase.append(0)

vecirc = QuantumCircuit(2*qubit_num+1)
vecirc.cx(0,3)
vec_elements.append(vecirc)
vec_coefficient.append(0.5*h)
vec_phase.append(0)

#3
vecirc = QuantumCircuit(2*qubit_num+1)
vecirc.cx(0,4)
vec_elements.append(vecirc)
vec_coefficient.append(0.5*h)
vec_phase.append(0)

#4
vecirc = QuantumCircuit(2*qubit_num+1)
vecirc.cx(0,1)
vecirc.cy(0,3)
vec_elements.append(vecirc)
vec_coefficient.append(0.125*gamma)
vec_phase.append(0)

#5
vecirc = QuantumCircuit(2*qubit_num+1)
vecirc.cx(0,3)
vecirc.cy(0,1)
vec_elements.append(vecirc)
vec_coefficient.append(0.125*gamma)
vec_phase.append(0)

#6
vecirc = QuantumCircuit(2*qubit_num+1)
vecirc.cx(0,1)
vecirc.cx(0,3)
vec_elements.append(vecirc)
vec_coefficient.append(0.125*gamma)
vec_phase.append(pi/2)

#7
vecirc = QuantumCircuit(2*qubit_num+1)
vecirc.cy(0,1)
vecirc.cy(0,3)
vec_elements.append(vecirc)
vec_coefficient.append(-0.125*gamma)
vec_phase.append(pi/2)

#9
vecirc = QuantumCircuit(2*qubit_num+1)
vecirc.cz(0,1)
vec_elements.append(vecirc)
vec_coefficient.append(-0.125*gamma)
vec_phase.append(pi/2)

#10
vecirc = QuantumCircuit(2*qubit_num+1)
vecirc.cz(0,3)
vec_elements.append(vecirc)
vec_coefficient.append(-0.125*gamma)
vec_phase.append(pi/2)

#4
vecirc = QuantumCircuit(2*qubit_num+1)
vecirc.cx(0,2)
vecirc.cy(0,4)
vec_elements.append(vecirc)
vec_coefficient.append(0.125*gamma)
vec_phase.append(0)

#5
vecirc = QuantumCircuit(2*qubit_num+1)
vecirc.cx(0,4)
vecirc.cy(0,2)
vec_elements.append(vecirc)
vec_coefficient.append(0.125*gamma)
vec_phase.append(0)

#6
vecirc = QuantumCircuit(2*qubit_num+1)
vecirc.cx(0,2)
vecirc.cx(0,4)
vec_elements.append(vecirc)
vec_coefficient.append(0.125*gamma)
vec_phase.append(pi/2)

#7
vecirc = QuantumCircuit(2*qubit_num+1)
vecirc.cy(0,2)
vecirc.cy(0,4)
vec_elements.append(vecirc)
vec_coefficient.append(-0.125*gamma)
vec_phase.append(pi/2)

#9
vecirc = QuantumCircuit(2*qubit_num+1)
vecirc.cz(0,2)
vec_elements.append(vecirc)
vec_coefficient.append(-0.125*gamma)
vec_phase.append(pi/2)

#10
vecirc = QuantumCircuit(2*qubit_num+1)
vecirc.cz(0,4)
vec_elements.append(vecirc)
vec_coefficient.append(-0.125*gamma)
vec_phase.append(pi/2)

vecirc = QuantumCircuit(2*qubit_num+1)
vec_elements.append(vecirc)
vec_coefficient.append(-0.5*gamma)
vec_phase.append(pi/2)



def generate_circuit_vec(index_diff, index_vec, ele):
    ancilla = QuantumRegister(1)
    phys = QuantumRegister(2*qubit_num)
    circ = QuantumCircuit(ancilla,phys)
    circ.h(0)
    circ.rz(vec_phase[index_vec],0)
    circ.x(0)
    for i in range(index_diff):
        circ.append(circuit_elements[i].to_instruction(),[phys[i] for i in range(2*qubit_num)])
    circ.barrier()
    circ.append(diff_elements[index_diff][ele].to_instruction(),list(range(2*qubit_num+1)))
    circ.barrier()
    circ.x(0)
    for i in range(index_diff,len(Parameters)):
        circ.append(circuit_elements[i].to_instruction(),[phys[i] for i in range(2*qubit_num)])
    circ.barrier()
    circ.append(vec_elements[index_vec].to_instruction(),list(range(2*qubit_num+1)))
    circ.barrier()
    circ.h(0)
    return circ
 
circuit_vec = [[] for i in range(len(Parameters))]
coeff_vec = [[] for i in range(len(Parameters))]

for i in range(len(Parameters)):
    circuit_vec[i] = []
    coeff_vec[i] = []
    for p,q in product(range(ele_num[i]),range(len(vec_elements))):
        circuit_vec[i].append(generate_circuit_vec(i,q,p))
        coeff_vec[i].append(diff_coeff[i][p]*vec_coefficient[q])
        
def evaluate_vec(circuit_vec, para, coeff_vec, shot_num):
    vec = np.ones(len(para),dtype = 'float')
    for i in range(len(Parameters)):
        value = 0
        for k, circ in enumerate(circuit_vec[i]):
            circuit = circ.bind_parameters({Parameters[m]: para[m] for m in range(len(Parameters))})
            circuit = transpile(circuit, basis_gates=basis_gates ,optimization_level=3)
            job = execute(circuit,backend = statevector_simulator)     
            result = job.result().get_statevector()
            p = 0 
            for m in range(2**(2*qubit_num)):
                  p+= abs(result[2*m])**2
            #print(i,k,(2*obtain_counts(counts,'0')/shot_num-1)*coeff_vec[i][k])
            value += (2*p-1)* coeff_vec[i][k]
        vec[i] = value
    return vec

vec = evaluate_vec(circuit_vec, para, coeff_vec, 2**16)

#dpara = np.linalg.solve(mat,vec)
#print(dpara)

def mybin(x,qubit_num):
    bstr = bin(x).replace('0b','')
    bstr = bstr.rjust(qubit_num,'0')
    return bstr

def measurement(para, shot_num):
    qc = QuantumCircuit(2)
    qc.ry(para[0],[0,1])
    job = execute(qc,backend = statevector_simulator)
    result = job.result().get_statevector()
    #print(result.data)
    rho = abs(result.data)
    rho = rho/sum(rho)
    #print(rho)
    
    p = np.zeros([2**qubit_num,2**qubit_num],dtype = complex)
    
    for i in range(2**qubit_num):
        qc = QuantumCircuit(2)
        bstr = mybin(i, qubit_num)
        #print(bstr)
        for k in range(qubit_num):
            if(bstr[k]=='1'):
                qc.x(k)
        qc.rx(para[1],[0,1])
        qc.ry(para[2],[0,1])
        qc.cz(0,1)
        #qc.rzz(para[3],0,1)
        qc.rx(para[3],[0,1])
        qc.ry(para[4],[0,1])
        qc.cz(0,1)
        
        job = execute(qc, backend = statevector_simulator)
        p[i] = job.result().get_statevector().data
        #print(p[i])
    mz = 0
    for i in range(2**qubit_num):
        mz += rho[i]*(abs(p[i][0])**2-abs(p[i][3])**2)
    return mz

#print(measurement(para,2**13))

# def zero_filter(mat):
#     for i,j in product(range(len(Parameters)),range(len(Parameters))):
#         #mat[i][j] = '%.6f' % mat[i][j]
#         if(abs(mat[i][j])<10**(-6)):
#             mat[i][j] = 0.
#     return mat


def diff_para(circuit_matrix, circuit_vec, coeff_matrix, coeff_vec, para, shot_num):
    mat = evaluate_matrix(circuit_matrix, para, coeff_matrix, shot_num)
    vec = evaluate_vec(circuit_vec, para, coeff_vec, shot_num)
    dpara = np.linalg.lstsq(mat, vec, rcond = 10**(-8))[0]

    return dpara, mat, vec

def TDVA_evolution(para, dt, steps, shot_num, circuit_matrix, circuit_vec, coeff_matrix, coeff_vec):
    mz = np.zeros(steps, dtype = float)
    t = np.zeros(steps, dtype = float)
    paras = np.zeros([steps,len(Parameters)], dtype = float)
    mats = np.zeros([steps,len(Parameters),len(Parameters)], dtype = float)
    vecs = np.zeros([steps,len(Parameters)], dtype = float)
    for i in range(steps):
        mz[i] = measurement(para,2**16)
        t[i] = i*dt
        paras[i] = para
        dpara, mat, vec = diff_para(circuit_matrix, circuit_vec, coeff_matrix, coeff_vec,para, shot_num)
        para += dpara*dt
        if(para[0]>pi):
            para[0] = pi
        if(para[0]<0):
            para[0] = 0
        mats[i] = mat
        vecs[i] = vec
        print(t[i],mz[i],para)
    return mz, t, paras, mats, vecs

mz, t, paras, mats, vecs = TDVA_evolution(para, 0.01, 500, 2**16, circuit_matrix, circuit_vec, coeff_matrix, coeff_vec)


font1 = {
    'family':'Times New Roman',
    'weight': 'normal',
    'size': 18,
    }

plt.figure(figsize = (12.0,9.0))
plt.tick_params(labelsize=20)

plt.plot(t, mz, label = 'ex state', color = 'black', linestyle = '-', marker = 'None')

plt.xlabel('γt', fontsize = 22)
plt.ylabel('value', fontsize = 22)
plt.legend(prop = font1, loc = 'center right')
figname = 'TFIMexact.jpg'
plt.savefig(figname)
plt.show()


dataname = 'TFIMexact.pkl'
f = open(dataname,'wb')
pickle.dump((mz, t, paras, mats, vecs),f)
f.close() 
