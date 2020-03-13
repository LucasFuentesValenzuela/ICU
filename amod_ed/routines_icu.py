#routines to support in the numerical tests for the routing algorithm
from amod_ed.helpers_icu import BPR, BPR_int
import numpy as np
import networkx as nx


#sign not necessary (think about what is actually D-1 and about the fact that it has a negative sign in front)
def update_costs(G):
    G=G.copy()
    UPPER_LIMIT = 10**22

    for e in G.edges:
        x = G[e[0]][e[1]]['f_m']+G[e[0]][e[1]]['f_r']
        phi = G[e[0]][e[1]]['phi']
        k = G[e[0]][e[1]]['k']
        G[e[0]][e[1]]['cost'] = BPR(phi, x, k)

        #if the capacity is too small (only happens for the edges to R)
        #then we can say that there is no flow there actually
        if k < 10**-5:  # we eliminate the edges with a too high cost from the equation
            #TODO: figure out whether this is a problem for the iterative alg?
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
            ri_k[e[0]] += G[e[0]][e[1]]['f_m']#Add to the origin node
            ri_k[e[1]] -= G[e[0]][e[1]]['f_m']#Substract from the destination node

    for n in ri_k.keys():
        if n!='R' and ri_k[n]>0: #node in deficit
            ri_k['R']+=ri_k[n]

    for n in G.nodes():
        if ri_prev==[]:
            G.nodes[n]["ri"]=ri_k[n]
        else:
            G.nodes[n]["ri"] = (1-beta) * ri_prev[n] + beta*ri_k[n]
        ri_k[n] = G.nodes[n]["ri"] 
    return ri_k, G


def update_OD(OD, ri_k, G, evolving_bounds=True):

    #update the OD pairs for rebalancers
    eps = 10**-6
    for n in ri_k.keys():
        if n != 'R' and not n.endswith('_p'):
            if ri_k[n] < -eps:  # you are in excess
                OD[(n, 'R')] = -ri_k[n]
            else:
                OD[(n, 'R')] = 0

    #TODO: improve the evolving bounds routine
    #currently, we keep lower bounds at zero
    if evolving_bounds:
        #TODO: parameters, to be given as arguments in the future
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
                # print("----")
                # print(crt_p_flow)
                # print(crt_Ulim)
                # print(rel_U_error)
                # UPPER BOUND
                if rel_U_error <= l1:  # x_k too close to U
                    new_Ulim = alpha*crt_Ulim
                elif rel_U_error >= l2:  # x_k too far from U
                    new_Ulim = alpha*crt_p_flow
                else:
                    new_Ulim = crt_Ulim
                # print(new_Ulim)
                OD[o, d] = new_Ulim
    return OD


def update_capacities(G, ri_k):
    eps = 10**-6
    for n in G.nodes():
        if not n.endswith('_p') and n != 'R': #rebalancing for those nodes means nothing
            if (n,'R') in G.edges:
                if ri_k[n] > eps:
                    G[n]['R']['k'] = ri_k[n]
                else:
                    G[n]['R']['k'] = eps
    return G

#we currently assign a rebalancing flow to y_k, but this is useless currently
#as with the new cost function we do not take it into account
def AoN(G, OD):
    #perform the All Or Nothing assignment
    y_k = init_y(G)
    eps = 10**-6

    # TODO: what is this alpha for? ANSWER: for the update below (in comments)
    # alpha = 1.5
    for (o, d) in OD.keys():
        U = OD[o, d]
        if U > eps:
            if d == 'R':
                flag = 'f_r'
            else:
                flag = 'f_m'
            path = nx.shortest_path(G, source=o, target=d, weight='cost')

            #TODO: figure out what the commented version here below actually meant
            #it seems to me like it was a fix to make sure we introduce the
            #"complement" of the demand
            #we have to see whether or not it was helping in the ICU case... ?
            if len(path) == 2 and d != 'R':
                # y_k[(o,d),'f_m']+=alpha*(U-G[dummy_nodes[d]][d]['f_m'])
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
    #The update factor quantifies to how much of the capacity we are allowed to go on any one rebalancing edge

    gamma = 2/(k+2) #initially only contained this

    if not update:
        return gamma

    #now we introduce a slight modification to avoid overshooting in the direction of too much rebalancing
    beta=[]#coefficients of updates
    for e in G_k.edges():
        if e[1]=='R':
            x_max=update_factor*G_k[e[0]][e[1]]['k']
            x_k_e=G_k[e[0]][e[1]]['f_m']+G_k[e[0]][e[1]]['f_r']
            y_k_e=y_k[(e[0],e[1]),'f_r']+y_k[(e[0],e[1]),'f_r']
            try:
                beta_=(x_max-x_k_e)/(y_k_e-x_k_e)
            except ZeroDivisionError:
                beta_=np.nan
            if beta_<0:
                beta.append(np.nan)
            else:
                beta.append(beta_)

    min_beta=np.nanmin(beta)
    if np.isnan(min_beta):
        return gamma    
    return np.minimum(gamma, min_beta) 


def update_flows(G, y_k, a_k, edge_list, solver='Outer'):
    #TODO: understand why the balance drops abruptly to near zero when reaching the last inner iteration
    #this is the case for ni = 50, 100, 1000, 5000,. ... and every time the value of the balance
    #drops to 10**-15, kind of magically, as if the last update was always perfect. 
    """
    We want to update the flows differently depending on using the ICU or the outer loop solver. 
    Why? Simply because in the ICU we do not use the rebalancing update for what we are doing. 
    """
    G=G.copy()
    if solver == 'Outer':
        flags = ['f_m', 'f_r']
    elif solver == 'ICU':
        flags = ['f_m']
    else:
        return

    # for i in range(len(edge_list)):
        # e = edge_list[i]
    for e in G.edges():
        for flag in flags:
            x_k_e = G[e[0]][e[1]][flag]  # retrieve the flow
            # retrieve the flow from the manual assignment
            y_k_e = y_k[(e[0], e[1]), flag]
           
            G[e[0]][e[1]][flag] = (1-a_k)*x_k_e + a_k * y_k_e
            # if e[0]=='0' and e[1]=='1':
                # print("in update flows")
                # print(x_k_e, y_k_e, flag, (1-a_k) * x_k_e, a_k * y_k_e, G[e[0]][e[1]][flag])
    return G