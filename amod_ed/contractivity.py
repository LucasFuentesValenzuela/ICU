import cvxpy as cp
from amod_ed.helpers_icu import BPR_int, BPR
import matplotlib.pyplot as plt
import numpy as np
import os

IMAGES_PATH = "/Users/lucasfuentes/ASL/Images"

def viz_costs(name, phi_p, phi_inv, k_p, k_inv, shift_inv):
    edges = ["Edge $(1,2)$", "Edge $(2,1)$"]
    #plot paths for both directions
    x = np.linspace(0, 10, 100)
    _, ax = plt.subplots(2,1,figsize = (10, 10))
    for i in range(2):
        cost = BPR(phi_p[i], x, k_p[i], beta = 4)
        inv_d = -BPR(phi_inv[i], x, k_inv[i], beta = 4) + shift_inv[i]
        ax[i].plot(x, cost, label = 'Cost')
        ax[i].plot(x, inv_d, label = 'Inverse Demand')
        ax[i].grid()
        ax[i].set_xlabel("Total Flow")
        ax[i].set_ylabel("Cost")
        ax[i].legend()
        ax[i].set_title(edges[i])
        ax[i].set_ylim([cost[0]-10, inv_d[0]+10])
    plt.savefig(os.path.join(IMAGES_PATH, name+".png"), transparent = True) 

def _construct_problem(phi_p, phi_inv, k_p, k_inv, shift_inv):
    f_p = cp.Variable(2)
    f_r = cp.Variable(2)
    r = cp.Parameter()
    constraints = [f_p>=0, f_r>=0, r == f_r[1]-f_r[0], f_p <= 10, f_r <=10]
    total_cost = 0
    for i in range(2):
        cost = BPR_int(phi_p[i], f_p[i]+f_r[i], k_p[i], beta = 4)
        inv_d = -BPR_int(phi_inv[i], f_p[i], k_inv[i], beta = 4) + shift_inv[i]*f_p[i]
        total_cost = total_cost + cost - inv_d
    obj = cp.Minimize(total_cost)
    prob = cp.Problem(obj, constraints)
    return f_p, f_r, r, prob

def sample_solutions(name ,phi_p, phi_inv, k_p, k_inv, shift_inv, nsamples=100, seed =0):

    f_p, _, r, prob = _construct_problem(phi_p, phi_inv, k_p, k_inv, shift_inv)

    np.random.seed(seed)
    samples = np.random.uniform(-10,10, nsamples)
    Tr = []
    for i in range(nsamples):
        r.value = samples[i]
        prob.solve(solver = cp.GUROBI)
        if prob.status!='optimal':
            print("iteration %d, status %s" %(i, prob.status))
        Tr.append(f_p.value[0] - f_p.value[1])

    dr= []
    dT = []
    for i in range(nsamples-1):
        dr.append(np.abs(samples[i]-samples[i+1]))
        dT.append(np.abs(Tr[i]-Tr[i+1]))
    plt.figure(figsize=(10,6))
    plt.scatter(dr, np.divide(dT, dr), s = 10, marker = 'o')
    # plt.plot(dr, np.ones(len(dr)), 'r')
    plt.grid()
    plt.ylabel('$\|Tr_1 - Tr_2\|/\|r_1-r_2\|$')
    plt.xlabel('$\|r_1-r_2\|$')
    plt.xlim([0, 20])
    plt.ylim([0,1.3])
    plt.savefig(os.path.join(IMAGES_PATH, name+".png"), transparent = True)    

    rat = np.divide(dT, dr)
    idx = np.where(rat>1)
    if idx[0].shape[0]>0:
        print("ratio larger than 1")
        for i in idx[0]: 
            print("     values of ri: ", samples[i], samples[i+1])
            print("     values of Tr: ", Tr[i], Tr[i+1])
    return dT, dr
