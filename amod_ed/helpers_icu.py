import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import cvxpy as cp


def BPR_int(phi, x, kappa, alpha=0.15, beta=4):
    return cp.multiply(phi,(x+cp.multiply(alpha/(beta+1),cp.power(x, (beta+1))/np.power(kappa, beta))))

#returns the value of BPR int, not just an expression as is the case above


def BPR_int_val(phi, x, kappa, alpha=0.15, beta=4):
    return phi*(x+alpha/(beta+1)*np.power(x, (beta+1))/np.power(kappa, beta))


def BPR(phi, x, kappa, alpha=0.15, beta=4):
    return phi*(1+alpha*(np.divide(x, kappa))**beta)

def BPR_cp(phi, x, kappa, alpha=0.15, beta=4):
    return phi*(1+alpha*(cp.power(np.divide(x, kappa), beta)))


def phi(l, t):
    return 36*l/t


def cost_per_edge(alpha, beta, phi_vec, flow_vec, kappa_vec, K_vec):
    c = BPR(alpha, beta, phi_vec, flow_vec, kappa_vec)-K_vec
    return c


#TODO: check this function!
def Value_Total_Cost(G):
    F_E = 0
    for e in G.edges():  # you know for sure exactly what edge it is for
        cost_edge=0
        x_k_e_m = G[e[0]][e[1]]['f_m']
        x_k_e_r = G[e[0]][e[1]]['f_r']

        #retrieve parameters to compute the BPR
        phi = G[e[0]][e[1]]['phi']
        k = G[e[0]][e[1]]['k']

        if k < 10**-5:  # you eliminate the edges that are considered non-usable
            G[e[0]][e[1]]['tot_cost']=np.nan
            continue
        # if e[1] == 'R':  # not including the cost of edges 1R and 2R might make sense, as we want to rebalance whatever happens
            # continue
            # pass #not sure we have to cancel it in the end

        # I am assuming there will be syntaxic problems there
        cost_edge += BPR_int_val(phi, x_k_e_m + x_k_e_r, k)

        #this has to be included because it is directly included in the definition of the cost function
        if G[e[0]][e[1]]['sign'] == (-1):  # we have a negative edge
            cost_edge -= (x_k_e_m + x_k_e_r)*G[e[0]][e[1]]['shift']  # INVERSE_DEMAND_SHIFT

        # not entirely sure this needs to be here
        if 'pot' in G.nodes[e[1]]:
            cost_edge += G.nodes[e[1]]['pot']*(x_k_e_m + x_k_e_r)

        F_E += cost_edge
        G[e[0]][e[1]]['tot_cost']=cost_edge
    return F_E, G







