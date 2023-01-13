# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 14:38:54 2022

@author: xinyu
"""

import numpy as np
from math import pi
import pickle
from matplotlib import pyplot as plt
from itertools import product
import time

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

def generate_hardware_efficient_ansatz(qubit_num, layer = 2):
    VecParameters = []
    for i,j in product(range(2*layer),range(qubit_num)):
        VecParameters.append(Parameter('θ'+str(qubit_num*i+j)))
    ele_num = np.ones(qubit_num*2*layer,dtype = int)

    circuit_elements =  [[]for i in range(qubit_num*2*layer)]
    
    for i,j in product(range(layer),range(qubit_num)):
        circuit_elements[qubit_num*(2*i)+j] = QuantumCircuit(qubit_num,name = 'rx'+str(j))
        circuit_elements[qubit_num*(2*i)+j].rx(VecParameters[qubit_num*(2*i)+j],[j])
        
        circuit_elements[qubit_num*(2*i+1)+j] = QuantumCircuit(qubit_num,name = 'ry'+str(j))
        circuit_elements[qubit_num*(2*i+1)+j].ry(VecParameters[qubit_num*(2*i+1)+j],[j])
    for i,k in product(range(layer),range(qubit_num-1)):
        circuit_elements[qubit_num*(2*i+1)+qubit_num-1].cz(k,k+1)
            
            
    diff_elements = [[] for i in range(qubit_num*2*layer)]
    diff_coeff = [[] for i in range(qubit_num*2*layer)]
    for i,j in product(range(2*layer),range(qubit_num)):
        diff_elements[qubit_num*i+j] = [[] for i in range(ele_num[qubit_num*i+j])]
        diff_coeff[qubit_num*i+j] = [[] for i in range(ele_num[qubit_num*i+j])]
        for k in range(ele_num[qubit_num*i+j]):
            if(i%2 == 0): 
                diff_elements[qubit_num*i+j][k] = QuantumCircuit(qubit_num+1,name = 'cx0'+str(j+1))
                diff_elements[qubit_num*i+j][k].cx(0,j+1)
            if(i%2 == 1):
                diff_elements[qubit_num*i+j][k] = QuantumCircuit(qubit_num+1,name = 'cy0'+str(j+1))
                diff_elements[qubit_num*i+j][k].cy(0,j+1)
            diff_coeff[qubit_num*i+j][k] = 0.5
    return VecParameters, ele_num, circuit_elements, diff_elements, diff_coeff

time_start = time.time()

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
R = 4

diagpara_num = R 
diagparas = np.array([1.,0.,0.,0.]) #diavecvecvecparas, the Mth value equals to 1 minus the sum of all previous M-1 value

layer = 2
vecpara_num = 2*qubit_num*layer
vecparas = np.zeros(vecpara_num,dtype = float)
para_num = diagpara_num+vecpara_num
para = np.hstack((diagparas,vecparas))  

VecParameters, ele_num, circuit_elements, diff_elements, diff_coeff = generate_hardware_efficient_ansatz(qubit_num,layer)

circuit_count_M = 0
circuit_count_V = 0


                        
def mybin(x,qubit_num):
    bstr = bin(x).replace('0b','')
    bstr = bstr.rjust(qubit_num,'0')
    return bstr

def evaluate_elements(statevector, qubit_num):
    p = 0
    for i in range(2**(qubit_num-1)):
        p += abs(statevector[2*i])**2
    return 2*p-1
    
def generate_circuit_matrix_Repp(index1, index2, ele1, ele2, p):
    bstr = mybin(p, qubit_num)
    ancilla = QuantumRegister(1)
    phys = QuantumRegister(qubit_num)
    circ = QuantumCircuit(ancilla, phys)
    
    for k in range(qubit_num):
        if(bstr[k]=='1'):
            circ.x(qubit_num-k)
    
    circ.h(0)
    circ.x(0)
    
    for i in range(index1):
        circ.append(circuit_elements[i].to_instruction(),[phys[i] for i in range(qubit_num)])
    circ.barrier()
    circ.append(diff_elements[index1][ele1].to_instruction(),list(range(qubit_num+1)))
    circ.barrier()
    circ.x(0)
    circ.barrier()
    
    if(index1<index2):
        for i in range(index1,index2):
            circ.append(circuit_elements[i].to_instruction(),[phys[i] for i in range(qubit_num)])
        circ.barrier()
        
    circ.append(diff_elements[index2][ele2].to_instruction(),list(range(qubit_num+1)))
    circ.barrier()
    circ.h(0)
    
    circ = transpile(circ, basis_gates = basis_gates,optimization_level=3)
    return circ

def generate_circuit_matrix_Repq(index, ele, p,q):
    bstrp = mybin(p, qubit_num)
    bstrq = mybin(q, qubit_num)
    
    ancilla = QuantumRegister(1)
    phys = QuantumRegister(qubit_num)
    circ = QuantumCircuit(ancilla, phys)
    
    circ.h(0)
    circ.rz(pi/2,0) ##this is an important part, and other forms of circuit could be used here
    
    for k in range(qubit_num):
        if(bstrq[k]=='1'):
            circ.cx(0,qubit_num-k)
    circ.x(0)
    
    for k in range(qubit_num):
        if(bstrp[k]=='1'):
            circ.cx(0,qubit_num-k)
    
    if(index>0):
        for i in range(index):
            circ.append(circuit_elements[i].to_instruction(),[phys[i] for i in range(qubit_num)])
    circ.barrier()
    circ.append(diff_elements[index][ele].to_instruction(),list(range(qubit_num+1)))
    circ.barrier()
    circ.x(0)    
    circ.h(0)
    
    circ = transpile(circ, basis_gates = basis_gates,optimization_level=3)
    return circ

def generate_circuit_matrix_Impq(index, ele, p,q):
    bstrp = mybin(p, qubit_num)
    bstrq = mybin(q, qubit_num)
    
    ancilla = QuantumRegister(1)
    phys = QuantumRegister(qubit_num)
    circ = QuantumCircuit(ancilla, phys)
    
    circ.h(0)
    circ.rz(pi/2,0)
    
    for k in range(qubit_num):
        if(bstrq[k]=='1'):
            circ.cx(0,qubit_num-k)
    circ.x(0)
    
    for k in range(qubit_num):
        if(bstrp[k]=='1'):
            circ.cx(0,qubit_num-k)
    
    if(index>0):
        for i in range(index):
            circ.append(circuit_elements[i].to_instruction(),[phys[i] for i in range(qubit_num)])
    circ.barrier()
    circ.append(diff_elements[index][ele].to_instruction(),list(range(qubit_num+1)))
    circ.barrier()
    circ.x(0)    
    circ.rx(pi/2,0)
    
    circ = transpile(circ, basis_gates = basis_gates,optimization_level=3)
    return circ

circuit_matrix_Repp = [[] for i in range(vecpara_num)]
circuit_matrix_Repq = [[] for i in range(vecpara_num)]
circuit_matrix_Impq = [[] for i in range(vecpara_num)]
coeff_matrix_Repp = [[] for i in range(vecpara_num)]
coeff_matrix_pq = [[] for i in range(vecpara_num)]

for i in range(len(VecParameters)):
    circuit_matrix_Repp[i] = [[[] for k in range(R)] for j in range(i+1)]
    circuit_matrix_Repq[i] = [[[] for k in range(R)] for j in range(R)]
    circuit_matrix_Impq[i] = [[[] for k in range(R)] for j in range(R)]
    coeff_matrix_Repp[i] = [[] for j in range(i+1)]
    
for i in range(len(VecParameters)):
    for j,k in product(range(i+1),range(R)):
        for p,q in product(range(ele_num[i]),range(ele_num[j])):
            circuit_matrix_Repp[i][j][k].append(generate_circuit_matrix_Repp(j,i,q,p,k))
            coeff_matrix_Repp[i][j].append(diff_coeff[j][q]*diff_coeff[i][p])
            circuit_count_M += 1

for i in range(len(VecParameters)):
    for j,k in product(range(R),range(R)):
        for p in range(ele_num[i]):
            circuit_matrix_Repq[i][j][k].append(generate_circuit_matrix_Repq(i,p,j,k))
            coeff_matrix_pq[i].append(diff_coeff[i][p])
            circuit_count_M += 1
    
for i in range(len(VecParameters)):
    for j,k in product(range(R),range(R)):
        for p in range(ele_num[i]):
            circuit_matrix_Impq[i][j][k].append(generate_circuit_matrix_Impq(i,p,j,k))
            circuit_count_M += 1
     
def evaluate_matrix_Repp(vecparas, shot_num):
    mat = np.ones([vecpara_num,vecpara_num,R],dtype = 'float')
    for i in range(vecpara_num):
        for j,k in product(range(i+1),range(R)):
            value = 0
            for l,circ in enumerate(circuit_matrix_Repp[i][j][k]):
                circuit = circ.bind_parameters({VecParameters[m]: vecparas[m] for m in range(i)})
                job = execute(circuit,backend = statevector_simulator)     
                result = job.result().get_statevector()
                value += evaluate_elements(result,qubit_num+1)*coeff_matrix_Repp[i][j][l]
            mat[i][j][k] = value
            mat[j][i][k] = value
    return mat

#mat_Repp = evaluate_matrix_Repp(vecparas, 0)

def evaluate_matrix_pq(vecparas, shot_num):
    matR = np.ones([vecpara_num,R,R], dtype = float)
    matI = np.ones([vecpara_num,R,R], dtype = float)
    for i,j,k in product(range(vecpara_num),range(R),range(R)):
        Rvalue = 0
        Ivalue = 0
        for l,circ in enumerate(circuit_matrix_Repq[i][j][k]):
            circuit = circ.bind_parameters({VecParameters[m]: vecparas[m] for m in range(i)})
            job = execute(circuit,backend = statevector_simulator)     
            result = job.result().get_statevector()
            Rvalue += evaluate_elements(result,qubit_num+1)*coeff_matrix_pq[i][l]
        matR[i][j][k] = Rvalue
        for l,circ in enumerate(circuit_matrix_Impq[i][j][k]):
            circuit = circ.bind_parameters({VecParameters[m]: vecparas[m] for m in range(i)})
            job = execute(circuit,backend = statevector_simulator)     
            result = job.result().get_statevector()
            Ivalue += evaluate_elements(result,qubit_num+1)*coeff_matrix_pq[i][l]
        matI[i][j][k] = Ivalue
    return matR, matI

#matR, matI = evaluate_matrix_pq(vecparas, 0)

def evaluate_matrix(diagparas, vecparas, shot_num):
    mat = np.zeros([R+vecpara_num,R+vecpara_num],dtype = float)
    matpp = evaluate_matrix_Repp(vecparas, shot_num)
    matRpq, matIpq = evaluate_matrix_pq(vecparas, shot_num)
    for i in range(R):
        mat[i][i] = 1
        for j in range(vecpara_num):
            value = 2*(diagparas[i]*matRpq[j][i][i])
            mat[i][j+R] = value
            mat[j+R][i] = value
    for i in range(vecpara_num):
        for j in range(i+1):
            vpp = matpp[i][j][:]
            vpp = vpp * (diagparas**2)
            value = sum(vpp)
            # if (i == j): 
            #     print(matpp[i][j][:])
            for p,q in product(range(R),range(R)):
                value += diagparas[p]*diagparas[q]*(matRpq[i][p][q]*matRpq[j][q][p]-matIpq[i][p][q]*matIpq[j][q][p]) 
            value = 2*value
            mat[i+R][j+R] = value
            mat[j+R][i+R] = value
    return mat
     
J = 1.
h = 1.
gamma = 0.1

vecirc_elements = []

#0
vecirc = QuantumCircuit(qubit_num+1)
vecirc_elements.append(vecirc)

#1
vecirc = QuantumCircuit(qubit_num+1)
vecirc.cz(0,1)
vecirc.cz(0,2)
vecirc_elements.append(vecirc)

#2
vecirc = QuantumCircuit(qubit_num+1)
vecirc.cx(0,1)
vecirc_elements.append(vecirc)

#3
vecirc = QuantumCircuit(qubit_num+1)
vecirc.cx(0,2)
vecirc_elements.append(vecirc)

#4
vecirc = QuantumCircuit(qubit_num+1)
vecirc.cy(0,1)
vecirc_elements.append(vecirc)

#5
vecirc = QuantumCircuit(qubit_num+1)
vecirc.cy(0,2)
vecirc_elements.append(vecirc)

#6
vecirc = QuantumCircuit(qubit_num+1)
vecirc.cz(0,1)
vecirc_elements.append(vecirc)

#7
vecirc = QuantumCircuit(qubit_num+1)
vecirc.cz(0,2)
vecirc_elements.append(vecirc)


vec_elements = []
vec_coefficient = []
vec1_coefficient = [] ##+1 or -1m due to the conjugate of sigma in B system
vecirc = [[]for i in range(2)]



#0
vec_elements.append([1,0])
vec_coefficient.append(J*1j)
vec1_coefficient.append(1)
#1
vec_elements.append([0,1])
vec_coefficient.append(-J*1j)
vec1_coefficient.append(1)
#2

vec_elements.append([2,0])
vec_coefficient.append(h*1j)
vec1_coefficient.append(1)

#3

vec_elements.append([3,0])
vec_coefficient.append(h*1j)
vec1_coefficient.append(1)

#4
vec_elements.append([0,2])
vec_coefficient.append(-h*1j)
vec1_coefficient.append(1)

#5
vec_elements.append([0,3])
vec_coefficient.append(-h*1j)
vec1_coefficient.append(1)

#6
vec_elements.append([0,0])
vec_coefficient.append(-1*gamma)
vec1_coefficient.append(1)

#7
vec_elements.append([2,2])
vec_coefficient.append(0.25*gamma)
vec1_coefficient.append(1)

#8
vec_elements.append([4,4])
vec_coefficient.append(-0.25*gamma)
vec1_coefficient.append(-1)

#9
vec_elements.append([2,4])
vec_coefficient.append(-0.25*gamma*1j)
vec1_coefficient.append(-1)

#10
vec_elements.append([4,2])
vec_coefficient.append(-0.25*gamma*1j)
vec1_coefficient.append(1)

#11
vec_elements.append([0,6])
vec_coefficient.append(-0.25*gamma)
vec1_coefficient.append(1)
#12

vec_elements.append([6,0])
vec_coefficient.append(-0.25*gamma)
vec1_coefficient.append(1)

#13
vec_elements.append([3,3])
vec_coefficient.append(0.25*gamma)
vec1_coefficient.append(1)

#14
vec_elements.append([5,5])
vec_coefficient.append(-0.25*gamma)
vec1_coefficient.append(-1)

#15
vec_elements.append([3,5])
vec_coefficient.append(-0.25*gamma*1j)
vec1_coefficient.append(-1)

#16
vec_elements.append([5,3])
vec_coefficient.append(-0.25*gamma*1j)
vec1_coefficient.append(1)

#17
vec_elements.append([0,7])
vec_coefficient.append(-0.25*gamma)
vec1_coefficient.append(1)

#18
vec_elements.append([7,0])
vec_coefficient.append(-0.25*gamma)
vec1_coefficient.append(1)

vec_coefficient = np.array(vec_coefficient)
vec1_coefficient = np.array(vec1_coefficient)
vec_coefficient = vec_coefficient*vec1_coefficient 

def generate_circuit_vec_Re1(index_vec, p, q):
    bstrp = mybin(p, qubit_num)
    bstrq = mybin(q, qubit_num)
    
    ancilla = QuantumRegister(1)
    phys = QuantumRegister(qubit_num)
    circ = QuantumCircuit(ancilla, phys)
    
    circ.h(0)
    
    for k in range(qubit_num):
        if(bstrq[k]=='1'):
            circ.cx(0,qubit_num-k)
    circ.x(0)
    
    for k in range(qubit_num):
        if(bstrp[k]=='1'):
            circ.cx(0,qubit_num-k)
     
    for i in range(len(VecParameters)):
        circ.append(circuit_elements[i].to_instruction(),[phys[i] for i in range(qubit_num)])
    circ.barrier()
    circ.append(vecirc_elements[index_vec].to_instruction(),list(range(qubit_num+1)))
    circ.x(0)
    circ.h(0)
    
    circ = transpile(circ, basis_gates = basis_gates,optimization_level=3)
    return circ

def generate_circuit_vec_Im1(index_vec, p, q):
    bstrp = mybin(p, qubit_num)
    bstrq = mybin(q, qubit_num)
    
    ancilla = QuantumRegister(1)
    phys = QuantumRegister(qubit_num)
    circ = QuantumCircuit(ancilla, phys)
    
    circ.h(0)
    
    for k in range(qubit_num):
        if(bstrq[k]=='1'):
            circ.cx(0,qubit_num-k)
    circ.x(0)
    
    for k in range(qubit_num):
        if(bstrp[k]=='1'):
            circ.cx(0,qubit_num-k)
     
    for i in range(len(VecParameters)):
        circ.append(circuit_elements[i].to_instruction(),[phys[i] for i in range(qubit_num)])
    circ.barrier()
    circ.append(vecirc_elements[index_vec].to_instruction(),list(range(qubit_num+1)))
    circ.x(0)
    circ.rx(pi/2,0)
    
    circ = transpile(circ, basis_gates = basis_gates,optimization_level=3)
    return circ

def generate_circuit_vec_Re2(index_vec, index_diff, ele_diff, p, q):
    bstrp = mybin(p, qubit_num)
    bstrq = mybin(q, qubit_num)
    
    ancilla = QuantumRegister(1)
    phys = QuantumRegister(qubit_num)
    circ = QuantumCircuit(ancilla, phys)
    
    circ.h(0)
    circ.rz(pi/2,0)
    
    for k in range(qubit_num):
        if(bstrq[k]=='1'):
            circ.cx(0,qubit_num-k)
    circ.x(0)
    
    for k in range(qubit_num):
        if(bstrp[k]=='1'):
            circ.cx(0,qubit_num-k)
     
    for i in range(index_diff):
        circ.append(circuit_elements[i].to_instruction(),[phys[i] for i in range(qubit_num)])
    circ.barrier()
    circ.append(diff_elements[index_diff][ele_diff].to_instruction(),list(range(qubit_num+1)))
    circ.x(0)
    for i in range(index_diff,len(VecParameters)):
         circ.append(circuit_elements[i].to_instruction(),[phys[i] for i in range(qubit_num)])
    circ.barrier()
    circ.append(vecirc_elements[index_vec].to_instruction(),list(range(qubit_num+1)))
    circ.barrier()
    
    circ.h(0)
    
    circ = transpile(circ, basis_gates = basis_gates,optimization_level=3)
    return circ
       
def generate_circuit_vec_Im2(index_vec, index_diff, ele_diff, p, q):
    bstrp = mybin(p, qubit_num)
    bstrq = mybin(q, qubit_num)
    
    ancilla = QuantumRegister(1)
    phys = QuantumRegister(qubit_num)
    circ = QuantumCircuit(ancilla, phys)
    
    circ.h(0)
    circ.rz(pi/2,0)
    
    for k in range(qubit_num):
        if(bstrq[k]=='1'):
            circ.cx(0,qubit_num-k)
    circ.x(0)
    
    for k in range(qubit_num):
        if(bstrp[k]=='1'):
            circ.cx(0,qubit_num-k)
     
    for i in range(index_diff):
        circ.append(circuit_elements[i].to_instruction(),[phys[i] for i in range(qubit_num)])
    circ.barrier()
    circ.append(diff_elements[index_diff][ele_diff].to_instruction(),list(range(qubit_num+1)))
    circ.x(0)
    for i in range(index_diff,len(VecParameters)):
         circ.append(circuit_elements[i].to_instruction(),[phys[i] for i in range(qubit_num)])
    circ.barrier()
    circ.append(vecirc_elements[index_vec].to_instruction(),list(range(qubit_num+1)))
    circ.barrier()
    
    circ.rx(pi/2,0)
    
    circ = transpile(circ, basis_gates = basis_gates,optimization_level=3)
    return circ

circuit_vec_Re1 = []
circuit_vec_Im1 = []
circuit_vec_Re2 = [[] for i in range(vecpara_num)]
circuit_vec_Im2 = [[] for i in range(vecpara_num)]
coeff_vec2 = [[] for i in range(vecpara_num)]

circuit_vec_Re1 = [[] for j in range(len(vecirc_elements))]
circuit_vec_Im1 = [[] for j in range(len(vecirc_elements))]
for i in range(vecpara_num):

    circuit_vec_Re2[i] = [[] for j in range(len(vecirc_elements))]
    circuit_vec_Im2[i] = [[] for j in range(len(vecirc_elements))]

for j in range(len(vecirc_elements)):
    circuit_vec_Re1[j] = [[[] for p in range(R)] for q in range(R)]
    circuit_vec_Im1[j] = [[[] for p in range(R)] for q in range(R)]
    
for i,j in product(range(vecpara_num), range(len(vecirc_elements))):
    circuit_vec_Re2[i][j] = [[[] for p in range(R)] for q in range(R)]
    circuit_vec_Im2[i][j] = [[[] for p in range(R)] for q in range(R)]

for j in range(len(vecirc_elements)):
    for p,q in product(range(R),range(R)):
        circuit_vec_Re1[j][p][q].append(generate_circuit_vec_Re1(j,p,q))
        circuit_vec_Im1[j][p][q].append(generate_circuit_vec_Im1(j,p,q))
        circuit_count_V += 2

for i,j in product(range(vecpara_num), range(len(vecirc_elements))):   
    for p,q,r in product(range(R),range(R),range(ele_num[i])):
        circuit_vec_Re2[i][j][p][q].append(generate_circuit_vec_Re2(j,i,r,p,q))      
        circuit_vec_Im2[i][j][p][q].append(generate_circuit_vec_Im2(j,i,r,p,q))   
        circuit_count_V += 2
        
for i in range(vecpara_num):
    for r in range(ele_num[i]):
        coeff_vec2[i].append(diff_coeff[i][r])
 
def evaluate_vec1(vecparas,shot_num):
    vecRe1 = np.zeros([len(vecirc_elements),R,R],dtype = float)
    vecIm1 = np.zeros([len(vecirc_elements),R,R],dtype = float)
    for j,p,q in product(range(len(vecirc_elements)),range(R),range(R)):
        circ = circuit_vec_Re1[j][p][q][0]
        circuit = circ.bind_parameters({VecParameters[m]: vecparas[m] for m in range(vecpara_num)})
        job = execute(circuit, backend = statevector_simulator)
        result = job.result().get_statevector()
        vecRe1[j][p][q] = evaluate_elements(result,qubit_num+1)
        
        circ = circuit_vec_Im1[j][p][q][0]
        circuit = circ.bind_parameters({VecParameters[m]: vecparas[m] for m in range(vecpara_num)})
        job = execute(circuit, backend = statevector_simulator)
        result = job.result().get_statevector()
        vecIm1[j][p][q] = evaluate_elements(result,qubit_num+1)
    
    return vecRe1, vecIm1

def evaluate_vec2(vecparas, shot_num):
    vecRe2 = np.zeros([vecpara_num,len(vecirc_elements),R,R],dtype = float)
    vecIm2 = np.zeros([vecpara_num,len(vecirc_elements),R,R],dtype = float)
    for i,j,p,q in product(range(vecpara_num),range(len(vecirc_elements)),range(R),range(R)):
        value = 0
        for l,circ in enumerate(circuit_vec_Re2[i][j][p][q]):
            circuit = circ.bind_parameters({VecParameters[m]: vecparas[m] for m in range(vecpara_num)})
            job = execute(circuit, backend = statevector_simulator)
            result = job.result().get_statevector()
            value += evaluate_elements(result,qubit_num+1)*coeff_vec2[i][l]
        vecRe2[i][j][p][q] = value
        
        value = 0
        for l,circ in enumerate(circuit_vec_Im2[i][j][p][q]):
            circuit = circ.bind_parameters({VecParameters[m]: vecparas[m] for m in range(vecpara_num)})
            job = execute(circuit, backend = statevector_simulator)
            result = job.result().get_statevector()
            value += evaluate_elements(result,qubit_num+1)*coeff_vec2[i][l]
        vecIm2[i][j][p][q] = value

    return vecRe2, vecIm2
 
def evaluate_vec(diagparas, vecparas, shot_num):
    vec = np.zeros(R+vecpara_num,dtype = float)
    vecRe1, vecIm1 = evaluate_vec1(vecparas,shot_num)
    vecRe2, vecIm2 = evaluate_vec2(vecparas,shot_num)
    vec1 = vecRe1+vecIm1*1j
    vec2 = vecRe2+vecIm2*1j
    for i in range(R):
        for j,p in product(range(len(vec_elements)),range(R)):
            value = vec1[vec_elements[j][0]][i][p]*vec1[vec_elements[j][1]][p][i]
            vec[i] += np.real(vec_coefficient[j]*diagparas[p]*value)
    for i in range(vecpara_num):
        for j,p,q in product(range(len(vec_elements)),range(R),range(R)):
            value = vec2[i][vec_elements[j][0]][p][q]*vec1[vec_elements[j][1]][q][p]+vec1[vec_elements[j][0]][p][q]*np.conj(vec2[i][vec_elements[j][1]][p][q])
            vec[i+R] += np.real(vec_coefficient[j]*diagparas[p]*diagparas[q]*value)
    return vec

time_end = time.time()
print('time cost for circuit preparation', time_end-time_start,'s')

time_start = time.time()
mat = evaluate_matrix(diagparas, vecparas, 0)
time_end = time.time()
print('number of circuit for matrix',circuit_count_M)
print('time cost for 1 evaluation of matrix', time_end-time_start,'s')
time_start = time.time()
vec = evaluate_vec(diagparas, vecparas, 0)      
time_end = time.time()
print('number of circuit for vector',circuit_count_V)
print('time cost for 1 evaluation of vector', time_end-time_start,'s')

# vec = evaluate_vec(circuit_vec, para, coeff_vec, 2**16)

#dpara = np.linalg.solve(mat,vec)
#print(dpara)

def measurement(diagparas, vecparas, shot_num):
 
    rho = diagparas/np.sum(diagparas)
    
    p = np.zeros([R,2**qubit_num],dtype = complex)
    
    for i in range(R):
        qc = QuantumCircuit(qubit_num)
        bstr = mybin(i, qubit_num)
        for k in range(qubit_num):
            if(bstr[k]=='1'):
                qc.x(qubit_num-k-1)
        for k in range(vecpara_num):
            qc.append(circuit_elements[k].to_instruction(),list(range(qubit_num)))
        
        qc = qc.bind_parameters({VecParameters[m]: vecparas[m] for m in range(vecpara_num)})
        
        job = execute(qc, backend = statevector_simulator)
        p[i] = job.result().get_statevector().data
        #print(p[i])
    mz = 0
    for i in range(R):
        mz += rho[i]*(abs(p[i][0])**2-abs(p[i][3])**2)
    
    return mz

# def zero_filter(mat):
#     for i,j in product(range(len(VecParameters)),range(len(VecParameters))):
#         #mat[i][j] = '%.6f' % mat[i][j]
#         if(abs(mat[i][j])<10**(-6)):
#             mat[i][j] = 0.
#     return mat


def diff_para(diagparas, vecparas, shot_num):
    mat = evaluate_matrix(diagparas, vecparas, shot_num)
    vec = evaluate_vec(diagparas, vecparas, shot_num)
    dpara = np.linalg.lstsq(mat, vec, rcond = 10**(-6))[0]
    return dpara, mat, vec

def TDVA_evolution(para, dt, steps, shot_num):
    mz = np.zeros(steps, dtype = float)
    t = np.zeros(steps, dtype = float)
    paras = np.zeros([steps,vecpara_num+R], dtype = float)
    mats = np.zeros([steps,vecpara_num+R,vecpara_num+R], dtype = float)
    vecs = np.zeros([steps,vecpara_num+R], dtype = float)
    for i in range(steps):
        mz[i] = measurement(para[:R],para[R:vecpara_num+R],2**16)
        t[i] = i*dt
        paras[i] = para
        dpara, mat, vec = diff_para(para[:R],para[R:vecpara_num+R], shot_num)
        para += dpara*dt
        mats[i] = mat
        vecs[i] = vec
        print(t[i],mz[i],para)
    return mz, t, paras, mats, vecs

mz, t, paras, mats, vecs = TDVA_evolution(para, 0.01, 500, 2**16)


font1 = {
    'family':'Times New Roman',
    'weight': 'normal',
    'size': 18,
    }

plt.figure(figsize = (12.0,9.0))
plt.tick_params(labelsize=20)
plt.plot(t, mz, label = 'average magnetization', color = 'black', linestyle = '-', marker = 'None')
plt.title('Average magnetization of TFIM ',fontsize=22)
plt.xlabel('t', fontsize = 22)
plt.ylabel('Average magnetization', fontsize = 22)
figname = 'TFIMexactver2.jpg'
plt.savefig(figname)
plt.show()


dataname = 'TFIMexactver2.pkl'
f = open(dataname,'wb')
pickle.dump((mz, t, paras, mats, vecs),f)
f.close() 
