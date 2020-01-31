import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import cvxpy as cp


def BPR_int(phi, x, kappa, alpha=0.15, beta=4):
    return phi*(x+alpha/(beta+1)*cp.power(x, (beta+1))/np.power(kappa, beta))

#returns the value of BPR int, not just an expression as is the case above


def BPR_int_val(phi, x, kappa, alpha=0.15, beta=4):
    return phi*(x+alpha/(beta+1)*np.power(x, (beta+1))/np.power(kappa, beta))


def BPR(phi, x, kappa, alpha=0.15, beta=4):
    return phi*(1+alpha*(np.divide(x, kappa))**beta)


def phi(l, t):
    return 36*l/t


def cost_per_edge(alpha, beta, phi_vec, flow_vec, kappa_vec, K_vec):
    c = BPR(alpha, beta, phi_vec, flow_vec, kappa_vec)-K_vec
    return c


def Value_Total_Cost(G):
    F_E = 0
    for e in G.edges():  # you know for sure exactly what edge it is for
        x_k_e_m = G[e[0]][e[1]]['f_m']
        x_k_e_r = G[e[0]][e[1]]['f_r']

        #retrieve parameters to compute the BPR
        phi = G[e[0]][e[1]]['phi']
        k = G[e[0]][e[1]]['k']

        if k < 10**-5:  # you eliminate the edges that are considered non-usable
            continue
        if e[1] == 'R':  # not including the cost of edges 1R and 2R might make sense, as we want to rebalance whatever happens
            continue

        # I am assuming there will be syntaxic problems there
        F_E += BPR_int_val(phi, x_k_e_m + x_k_e_r, k)

        #this has to be included because it is directly included in the definition of the cost function
        if G[e[0]][e[1]]['sign'] == (-1):  # we have a negative edge
            F_E -= (x_k_e_m + x_k_e_r)*80  # INVERSE_DEMAND_SHIFT

        # not entirely sure this needs to be here
        # if 'pot' in G.nodes[e[1]]:
        #     F_E+=G.nodes[e[1]]['pot']*flow_tmp

    return F_E







