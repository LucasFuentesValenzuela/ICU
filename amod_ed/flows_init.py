import cvxpy as cp
import numpy as np
import networkx as nx
from amod_ed.routines_icu import update_costs


def initialize_flows(G, G_prev, ri_k, OD):
    """
    Initialize the flows for graph G based on the flows of the last inner
    iteration. 
    """

    # init the flows to zero
    for e in G.edges():
        for flag in ['f_r', 'f_m']:
            G[e[0]][e[1]][flag] = 0

    if G_prev == None:
        # in this case, assign by shortest path as before
        G = init_flows_shortestPath(G, OD, G_prev=None)

    else:
        # keep the flows from the previous network (passengers)
        G = initialize_passengers(G, G_prev)
        # determine the best start for the rebalancers
        G = initialize_rebalancers(G, G_prev, ri_k)

    return G


def initialize_rebalancers(G, G_prev, ri_k):
    """
    Rebalancers are initialized by finding the closest feasible point
    for the rebalancing flow f_r (the passenger flow is feasible as OD don't change). 
    """

    r_i, nodes_list = _get_rebalancing_vector(ri_k)
    f_r, edge_list = _get_rebalancing_flows(G, G_prev)
    # here we have G as it is updated with capacities
    A = _get_matrix_A(G, r_i, edge_list, nodes_list)
    B = _get_matrix_B(r_i, edge_list, nodes_list, A)
    f_init_r = _compute_nearest_feasible(f_r, r_i, A, B, norm=2)
    G = introduce_rebalancers(G, f_init_r, edge_list)

    return G


def initialize_passengers(G, G_prev):
    """
    Initializing passengers amounts to just taking the previous values of the
    passenger flows across the graph. 
    """
    # We keep the flows of the passengers.
    for e in G.edges():
        G[e[0]][e[1]]['f_m'] = G_prev[e[0]][e[1]]['f_m']
    return G


def _get_rebalancing_vector(ri_k):
    """
    Extract the rebalancing vector from the graph to compute the initialization. 

    Out
    ---
    vector of dimension n_nodes where each entry is the r_i
    Reminder: 
        r_i < 0 if node i in excess of rebalancers
        r_i > 0 if node i in deficit of rebalancers
    """
    r_i=[]
    nodes_list_all = list(ri_k.keys())
    nodes_list=[]
    for i in range(len(nodes_list_all)):
        n=nodes_list_all[i]

        if n.endswith('_p'): #This is a "dummy node"
            continue
        #we only keep the nodes that are actually in the graph + the rebalancing node

        nodes_list.append(n)
        r_i.append(ri_k[n])

    return np.array(r_i), nodes_list #now nodes_list is only the nodes that we keep


def _get_rebalancing_flows(G, G_prev):
    """
    Extract the rebalancing flows to compute the initialization
    """
    eps=10**-5
    edge_list_all = list(G_prev.edges())
    f_r = [] 
    edge_list = []
    for i in range(len(edge_list_all)):
        e = edge_list_all[i]

        if e[1].endswith('_p'): #this is a dummy edge
            continue
        if G[e[0]][e[1]]['k']<eps and e[1]=='R':#this edge is "deactivated" 
            continue

        edge_list.append(e)
        f_r.append(G_prev[e[0]][e[1]]['f_r'])
    return np.array(f_r), edge_list

def _get_matrix_B(r_i, edge_list, nodes_list, A):
    """
    This matrix enforces rebalancing flow conservation at each node with r_i >0. 

    Another way to do it (easier) is just to say that f_r (n to R) = r_i[n]
    """
    eps = 10**-5
    n_edges = len(edge_list)
    n_nodes = len(nodes_list)
    B = np.zeros((0, n_edges)) #create empty array

    for i in range(n_nodes):
        if r_i[i] > eps: 
            idx_e = edge_list.index((nodes_list[i], 'R'))
            new_row = np.zeros((1, n_edges))
            new_row[0, idx_e] = -1
            new_row = new_row + A[i,:]
            B = np.concatenate((B, new_row), axis = 0)
    
    if B.shape[0]==0:
        return np.zeros((1, n_edges))

    return B

def _get_matrix_A(G, r_i, edge_list, nodes_list):
    """
    Extract the out matrix from the graph to compute initialization. 

    Out
    ---
    Matrix n_nodes x n_edges, where entry ij is 1 only if i is origin of edge j
    """
    eps = 10**-5
    n_edges = len(edge_list)
    n_nodes = len(nodes_list)
    A_out = np.zeros((n_nodes, n_edges))
    A_in = np.zeros((n_nodes, n_edges))

    
    for i in range(len(nodes_list)):
        for j in range(len(edge_list)):
            n = nodes_list[i]
            e = edge_list[j]
            
            # equivalent to saying that the node is in excess of rebalancers, ie ri<0
            if G[e[0]][e[1]]['k'] < eps and e[1]=='R':  
                continue

            if n == e[0] and e[1]!='R':  # n is origin and is a graph edge
                A_out[i, j] = 1
            elif n == e[1] and n!='R':  # n is destination
                A_in[i, j] = 1

    A = A_in - A_out
    #currently the A matrix should have one row of zeros -- TODO: discard that

    # for i in range(len(nodes_list)):
    #     if r_i[i] > eps: 
    #         A[idx_R, :] = A[idx_R, :] - A[i,:]
    return A


def _compute_nearest_feasible(f_r, r_i, A, B, norm=2):
    """
    We compute the nearest feasible point for the rebalancing flow. 
    """
    f = cp.Variable(f_r.shape[0])
    constraints = [A*f == r_i, B*f == 0, f >= 0]
    obj = cp.Minimize(cp.norm(f-f_r, norm))
    prob = cp.Problem(obj, constraints)
    _ = prob.solve(solver = cp.ECOS)
    print("Initialization problem status: ", prob.status)
    f_init_r = f.value
    return f_init_r


def introduce_rebalancers(G_k, f_init_r, edge_list):
    """
    Introduce the new values of rebalancing flow computed by the initialization 
    into the graph
    """
    for i in range(len(edge_list)):
        e = edge_list[i]
        G_k[e[0]][e[1]]['f_r'] = f_init_r[i]
    return G_k


def init_flows_shortestPath(G, OD, G_prev=None):
    """
    Initialize the flow with a feasible solution.
    We keep the progress made in the last iteration. 
    """
    # TODO: check if this version works on ICU, too (minor change, see version in FW_ICU.py)

    # reinitialize the flows to zero
    for e in G.edges():
        for flag in ['f_r', 'f_m']:
            G[e[0]][e[1]][flag] = 0

    if not G_prev == None:
        """
        We keep the progress made by the previous iteration for the passenger flows. 
        Executed if we pass a previous G_prev as argument
        """
        # We keep the flows of the passengers.
        for e in G.edges():
            for flag in ['f_m']:
                G[e[0]][e[1]][flag] = G_prev[e[0]][e[1]][flag]
        # we update the costs so that the rebalancers are appropriately assigned there.
        G = update_costs(G)
        # We only update the flows of the rebalancers.
        for (o, d) in OD.keys():
            if d == 'R':
                path = nx.shortest_path(G, source=o, target=d, weight='cost')
                for i in range(len(path)-1):
                    G[path[i]][path[i+1]]['f_r'] += OD[o, d]

    # no G_prev given, just initialize as planned
    else:
        for (o, d) in OD.keys():
            path = nx.shortest_path(G, source=o, target=d, weight='cost')
            for i in range(len(path)-1):
                if d == 'R':
                    G[path[i]][path[i+1]]['f_r'] += OD[o, d]
                else:
                    G[path[i]][path[i+1]]['f_m'] += OD[o, d]
    return G
