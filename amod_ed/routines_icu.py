# routines to support in the numerical tests for the routing algorithm
from amod_ed.helpers_icu import BPR, BPR_int
import numpy as np
import networkx as nx
import cvxpy as cp


def update_costs(G):
    G = G.copy()
    UPPER_LIMIT = 10**22

    for e in G.edges:
        x = G[e[0]][e[1]]['f_m']+G[e[0]][e[1]]['f_r']
        phi = G[e[0]][e[1]]['phi']
        k = G[e[0]][e[1]]['k']
        G[e[0]][e[1]]['cost'] = BPR(phi, x, k)

        # if the capacity is too small (only happens for the edges to R)
        # then we can say that there is no flow there actually
        if k < 10**-5:
            G[e[0]][e[1]]['cost'] = UPPER_LIMIT
            G[e[0]][e[1]]['f_r'] = 0
            continue

        if G[e[0]][e[1]]['sign'] == (-1):  # we have a negative edge
            G[e[0]][e[1]]['cost'] -= G[e[0]][e[1]]['shift']

        if 'pot' in G.nodes[e[1]]:
            G[e[0]][e[1]]['cost'] += G.nodes[e[1]]['pot']
    return G


def estimate_ri_k(G, ri_smoothing, a_k, ri_prev):
    """
    Determine whether each node is in excess or deficit of rebalancers.
    """
    ri_k = dict()

    if ri_smoothing:
        beta = a_k
    else:
        beta = 1  # no smoothing

    for n in G.nodes():
        ri_k[n] = 0
        # ri_prev[n] = G.nodes[n]["ri"]
    for e in G.edges():
        if not e[1].endswith('_p') and e[1] != 'R':
            ri_k[e[0]] += G[e[0]][e[1]]['f_m']  # Add to the origin node
            # Substract from the destination node
            ri_k[e[1]] -= G[e[0]][e[1]]['f_m']

    for n in G.nodes():
        if ri_prev == []:
            G.nodes[n]["ri"] = ri_k[n]
        else:
            G.nodes[n]["ri"] = (1-beta) * ri_prev[n] + beta*ri_k[n]
        ri_k[n] = G.nodes[n]["ri"]

    return ri_k, G


def update_OD(OD, ri_k, G, evolving_bounds=False):

    # update the OD pairs for rebalancers
    eps = 10**-6
    for n in ri_k.keys():
        if n != 'R' and not n.endswith('_p'):
            if ri_k[n] < -eps:  # you are in excess
                OD[(n, 'R')] = -ri_k[n]
            else:
                OD[(n, 'R')] = 0

    # TODO: improve the evolving bounds routine
    # currently, we keep lower bounds at zero
    if evolving_bounds:
        # TODO: parameters, to be given as arguments in the future
        l1 = 10**-3
        l2 = 0.1  # such a high value currently disables it
        alpha = 1.2  # I think a moving value on that would be better. And same for the l1/l2. should be based o2
        for (o, d) in OD.keys():
            if d.endswith('_p'):  # only the nodes that are the dummy nodes
                crt_Ulim = OD[o, d]  # currently, we treat only the upper limit
                crt_p_flow = get_total_flow_to_dummy_node(G, d)
                if crt_p_flow == 0:
                    continue
                rel_U_error = abs(crt_p_flow-crt_Ulim)/crt_p_flow
                if rel_U_error <= l1:  # x_k too close to U
                    new_Ulim = alpha*crt_Ulim
                elif rel_U_error >= l2:  # x_k too far from U
                    new_Ulim = alpha*crt_p_flow
                else:
                    new_Ulim = crt_Ulim
                OD[o, d] = new_Ulim
    return OD


def update_capacities(G, ri_k):
    eps = 10**-6
    for n in G.nodes():
        # rebalancing for those nodes means nothing
        if not n.endswith('_p') and n != 'R':
            if (n, 'R') in G.edges:
                if ri_k[n] > eps:
                    G[n]['R']['k'] = ri_k[n]
                else:
                    G[n]['R']['k'] = eps
    return G


def AoN(G, OD):
    # perform the All Or Nothing assignment
    y_k = init_y(G)
    eps = 10**-6

    for (o, d) in OD.keys():
        U = OD[o, d]
        if U > eps:
            if d == 'R':
                flag = 'f_r'
            else:
                flag = 'f_m'
            path = nx.shortest_path(G, source=o, target=d, weight='cost')

            if len(path) == 2 and d != 'R':
                y_k[(o, d), 'f_m'] += U
            else:
                for i in range(len(path)-1):
                    y_k[(path[i], path[i+1]), flag] += U
    return y_k


def init_y(G):
    y_k = dict()
    for e in G.edges():
        y_k[e, 'f_r'] = 0
        y_k[e, 'f_m'] = 0
    return y_k


def get_total_flow_to_dummy_node(G, d):
    """
    Computes the total flow going into a dummy node. 
    This is useful when updating the bound, as the bound is a measure of the max flow that can go into
    a dummy node. 
    """
    flow = 0

    for e in G.edges():
        o_ = d.split('_')[0]
        if e[1] == d and e[0] != o_:
            flow += G[e[0]][d]['f_m']
    return flow


########################################
#
# Consolidated version
#
#######################################


def fixed_step(k, y_k, G_k, update=True, update_factor=1.1):
    """
    Simplest version of the fixed step. I think I made a big mistake in implementing the previous ones (see main files). 
    """
    # The update factor quantifies to how much of the capacity we are allowed to go on any one rebalancing edge

    gamma = 2/(k+2)  # initially only contained this

    if not update:
        return gamma

    # now we introduce a slight modification to avoid overshooting in the direction of too much rebalancing
    beta = []  # coefficients of updates
    for e in G_k.edges():
        if e[1] == 'R':
            x_max = update_factor*G_k[e[0]][e[1]]['k']
            x_k_e = G_k[e[0]][e[1]]['f_m']+G_k[e[0]][e[1]]['f_r']
            y_k_e = y_k[(e[0], e[1]), 'f_r']+y_k[(e[0], e[1]), 'f_r']
            try:
                beta_ = (x_max-x_k_e)/(y_k_e-x_k_e)
            except ZeroDivisionError:
                beta_ = np.nan
            if beta_ < 0:
                beta.append(np.nan)
            else:
                beta.append(beta_)

    min_beta = np.nanmin(beta)
    if np.isnan(min_beta):
        return gamma
    return np.minimum(gamma, min_beta)


def update_flows(G, y_k, a_k, edge_list, solver='Outer'):
    """
    We want to update the flows differently depending on using the ICU or the outer loop solver. 
    Why? Simply because in the ICU we do not use the rebalancing update for what we are doing. 
    """
    G = G.copy()
    if solver == 'Outer':
        flags = ['f_m', 'f_r']
    elif solver == 'ICU':
        flags = ['f_m']
    else:
        return

    for e in G.edges():
        for flag in flags:
            x_k_e = G[e[0]][e[1]][flag]  # retrieve the flow
            # retrieve the flow from the manual assignment
            y_k_e = y_k[(e[0], e[1]), flag]
            G[e[0]][e[1]][flag] = (1-a_k)*x_k_e + a_k * y_k_e
    return G


##################################################
#
# Line Search
#
##################################################




def line_search(G, y_k, edge_list):

    a_k = cp.Variable()

    constraints = [a_k >= 0, a_k <= 1]
    obj = Total_Cost_line_search(G, y_k, a_k, edge_list)
    # print(obj)

    prob = cp.Problem(cp.Minimize(obj), constraints) 
    prob.solve(solver=cp.ECOS, verbose=False)
    print('solver ECOS:', a_k.value)
    print("Status of line search problem : ", prob.status)

    # prob.solve(solver = cp.GUROBI, verbose=False)
    # print('solver GUROBI', a_k.value)
    # print("Status of line search problem : ", prob.status)

    # prob = cp.Problem(cp.Minimize(obj), constraints)
    # prob.solve(solver = cp.CVXOPT, verbose=False)
    # print('solver CVXOPT', a_k.value)
    # print("Status of line search problem : ", prob.status)

    return a_k.value

#TODO: rewrite everything as there is a lot of uncertainty around
#what should and should not be included
def Total_Cost_line_search(G, y_k, a_k, edge_list):

    F_E = 0

    for e in edge_list:  # you know for sure exactly what edge it is for

        cost_edge = 0
        x_k_e_m = G[e[0]][e[1]]['f_m']
        x_k_e_r = G[e[0]][e[1]]['f_r']
        flow_x = x_k_e_m + x_k_e_r

        y_k_e_m = y_k[(e[0], e[1]), 'f_m']
        y_k_e_r = y_k[(e[0], e[1]), 'f_r']
        flow_y = y_k_e_m + y_k_e_r

        delta = flow_y - flow_x

        #TODO: is this useful? 
        # guard against instabilities:
        # if np.abs(flow_x) < 10**-5:
        #     flow_x = 0
        # if np.abs(delta) < 10**-5:
        #     delta = 0

        flow_tmp = flow_x + a_k*delta
        print(e, flow_tmp)
        # retrieve parameters to compute the BPR
        phi = G[e[0]][e[1]]['phi']
        k = G[e[0]][e[1]]['k']

        #TODO: what exactly do you need here? 
        if k < 10**-5:  # you eliminate the edges that are considered non-usable
            continue
        # if e[1] == 'R':
            # continue

        cost_edge += BPR_int(phi, flow_tmp, k, beta=4)
        if G[e[0]][e[1]]['sign'] == (-1):  # we have a negative edge
            cost_edge -= (flow_tmp)*G[e[0]][e[1]]['shift']
        if 'pot' in G.nodes[e[1]]:
            cost_edge += G.nodes[e[1]]['pot']*(flow_tmp)

        F_E += cost_edge
    return F_E


###############################
#
# Flow conservation at nodes
#

def check_flow_cons(G, OD): 
    """ 
    Loop through the nodes and compute the total incoming and outgoing flow. 
    It should cancel at every origin/destination node if the network is properly rebalanced. 
    For other nodes (R and n_p) it should be equal to the total flow assigned to them
    in OD.
    """
    net_flow=dict()
    for n in G.nodes():
       net_flow[n]=0 
    for e in G.edges():
        o=e[0]
        d=e[1]
        if d.endswith('_p') or d=='R':#we are on a dummy edge on the graph
            continue
        net_flow[o]-=G[o][d]['f_m']+G[o][d]['f_r']
        net_flow[d]+=G[o][d]['f_m']+G[o][d]['f_r']
    
    #the set of nodes for which the flows should indeed be balanced. 
    l=[]
    for (o,d) in OD.keys():
        l.append(o)
        if d.endswith('_p'):
            d=d.split('_')[0]
        if d=='R':
            continue
        l.append(d)
    l=np.unique(np.array(l))
    return net_flow, l

def check_flow_cons_at_OD_nodes(G,OD):
    """
    Compute the flow balance only at nodes that belong to OD. 
    """
    net_flow,l=check_flow_cons(G,OD)
    balance=[]
    for n in l:
        balance.append(net_flow[n])
    return np.array(balance)