import cvxpy as cp
from amod_ed.helpers_icu import BPR_int, BPR, BPR_cp
import matplotlib.pyplot as plt
import numpy as np
import os
import networkx as nx


"""
A more complex and complete version of the contractivity framework
"""
IMAGES_PATH = "/Users/lucasfuentes/ASL/Images"


def viz_costs(edges, inv_edges, name = '', save = False):
    """
    Visualize the cost functions. 

    Parameters
    ----------
    edges: pd.DataFrame
        dataframe containing the edges
    inv_edges: pd.DataFrame 
        dataframe containing the inverse edges
    save: bool
        whether or not to save
    """
    G = nx.DiGraph()

    for i in range(len(edges)):
        G.add_edge(edges.loc[i, 'origin'], edges.loc[i, 'destination'])

    #plot paths for both directions
    x = np.linspace(0, 10, 100)
    _, ax = plt.subplots(len(inv_edges)+1,1,figsize = (10, 5*(len(inv_edges)+1)))

    #loop over both edges
    for i in range(len(inv_edges)):
        """
        Loop over od pairs
        """
        o = inv_edges.loc[i, 'origin']
        d = inv_edges.loc[i, 'destination']
        phi_inv = inv_edges.loc[i, 'phi']
        k_inv = inv_edges.loc[i, 'k']
        shift_inv = inv_edges.loc[i, 'shift']
        for p in nx.all_simple_paths(G, o, d):
            cost = 0
            for j in range(len(p)-1):
                eo = p[j]
                ed = p[j+1]
                e_ind = edges.index[np.logical_and(
                    edges['origin'] == eo, edges['destination'] == ed
                )][0]
                phi_p = edges.loc[e_ind, 'phi']
                k_p = edges.loc[e_ind, 'k']
                cost += BPR(phi_p, x, k_p, beta = 4)


            ax[i].plot(x, cost, label = str(p))
        inv_d = -BPR(phi_inv, x, k_inv, beta = 4) + shift_inv
        ax[i].plot(x, inv_d, label = 'Inverse Demand')
        ax[i].grid()
        ax[i].set_xlabel("Total Flow")
        ax[i].set_ylabel("Cost")
        ax[i].legend()
        ax[i].set_title('OD (%d, %d)'%(o,d))
        ax[i].set_ylim([cost[0]-10, inv_d[0]+10])
    if save:
        plt.savefig(os.path.join(IMAGES_PATH, name+".png"), transparent = True) 

def _construct_problem(edges, inv_edges):
    """
    Build the optimization problem based on the shape of the cost functions. 

    Parameters
    ----------
    edges: pd.DataFrame
        dataframe containing the edges
    inv_edges: pd.DataFrame 
        dataframe containing the inverse edges
    """
    G = nx.DiGraph()
    for i in range(len(edges)):
        G.add_edge(
            edges.loc[i, 'origin'], edges.loc[i, 'destination']
            )
    nodes = list(G.nodes)
    edge_list = list(G.edges())
    inv_edge_list = []
    for j in range(len(inv_edges)):
        inv_edge_list.append(
            (inv_edges.loc[j, 'origin'], inv_edges.loc[j, 'destination']))


    #construct a map for the coordinates
    map_comps = dict()
    map_edges = dict()
    for i in range(len(edge_list)):
        n_e = []
        e = edge_list[i]
        for j in range(len(inv_edge_list)):
            n = i*len(inv_edges)+j
            e_inv = inv_edge_list[j]
            map_comps[n] = [e, e_inv]
            n_e.append(n)
        map_edges[e] = n_e


    
    #Define the optimization variables
    f_p = cp.Variable(len(edges)*len(inv_edges)) 
    f_r = cp.Variable(len(edges))
    r = cp.Parameter(len(nodes))#in our setting, the rebalacing is a parameter
    d_var = cp.Variable(len(inv_edges))

    constraints = [
        f_p>=0, 
        f_r>=0, 
        f_p <= 10, 
        f_r <=10
        ]
    
    #constraint on r
    for i in range(len(nodes)):
        n = nodes[i]
        c_vec =[]
        for j in range(len(edges)):
            e = edges.loc[j, 'origin'], edges.loc[j, 'destination'] 
            if e[0] == n:
                c_vec.append(1)
            elif e[1] == n:
                c_vec.append(-1)
            else:
                c_vec.append(0)
        constraints.append(
            np.array(c_vec)*f_r == r[i]
        )
        
    # constraint on the flow 
    for i in range(len(inv_edge_list)):
        (o,d) = inv_edge_list[i]
        idx_end = []
        for p in nx.all_simple_paths(G, o, d):
            idx_list = []
            for j in range(len(p)-1):
                e = (p[j], p[j+1])
                for k in  map_edges[e]:
                    if map_comps[k] == [e, (o,d)]:
                        idx_list.append(k)
            for j in range(len(idx_list)-1):
                ii = idx_list[j]
                ii_next = idx_list[j+1]
                constraints.append(
                    f_p[ii] == f_p[ii_next]
                )
            idx_end.append(idx_list[-1])
        
        con_d = 0
        for ii in idx_end:
            con_d += f_p[ii]
        constraints.append(con_d == d_var[i])


    total_cost = 0
    #iterate over edges to build the total cost
    #build the edge costs
    costs_dict = dict()
    for i in range(len(edges)):
        #For each direction, there are two terms to the cost
        #1. the actual cost to the passenger
        #2. The inverse demand cost
        e = edges.loc[i, 'origin'], edges.loc[i, 'destination'] 
        phi_p = edges.loc[i, 'phi']
        k_p = edges.loc[i,'k']

        #build the total cost
        idx_e = map_edges[e]
        f_p_e = 0
        for ii in idx_e:
            f_p_e += f_p[ii]
        cost = BPR_int(phi_p, f_p_e+f_r[i], k_p, beta = 4)
        #The total cost is the sum of all those terms
        total_cost = total_cost + cost

        costs_dict[e] = BPR_cp(phi_p,f_p_e+f_r[i], k_p, beta = 4)

    inv_demand_dict = dict()
    for i in range(len(inv_edges)):
        ie = inv_edge_list[i]
        k_inv = inv_edges.loc[i, 'k']
        phi_inv = inv_edges.loc[i, 'phi']
        shift = inv_edges.loc[i, 'shift']
        inv_d = -BPR_int(phi_inv, d_var[i], k_inv, beta = 4) + shift*d_var[i]
        total_cost = total_cost - inv_d
        inv_demand_dict[ie] = - BPR_cp(phi_inv, d_var[i], k_inv, beta = 4) + shift

    #The objective is to minimize the total cost
    obj = cp.Minimize(total_cost)
    prob = cp.Problem(obj, constraints)


    """
    Returns
    -------
    f_p: cvxpy.Variable
        The flow for each commodity on each edge
    f_r: cvxpy.Variable
        The rebalancing flow on each edge
    r: cvxpy.Parameter
        The rebalancing guess for each node
    d_var: cvxpy.Variable
        The demand for each each
    prob: cvxpy.Problem
        The optimization problem
    map_comps: dict
        A map linking components of f_p to the edges and inv edges
    map_edges: dict
        A map linking edges to components of f_p
    costs_dict: dict
        Dict containing the cost for each edge
    inv_demand_dict: dict
        The inverse demand cost for each od pair
    G: nx.DiGraph
        Graph representing the network
    nodes: list
        list of nodes
    """
    return f_p, f_r, r, d_var, prob, map_comps, map_edges, costs_dict, inv_demand_dict, G, nodes





def run_algorithm(edges, inv_edges, nsolutions = 5, seed =0, max_iter = 50):
    """
    Run the algorithm starting from a series of randomly
    sampled initial points

    Parameters
    ----------
    nsolutions: int
        number of initial points to be drawn, from which to iterate
    seed: int
        seed for the random number generator
    max_iter: int
        maximum number of iterations after which to stop the algorithm

    Returns
    -------
    r_tot: list of lists
        each list contains the different values of r over the iterations
    """ 

    #build the optimization problem
    f_p, f_r, r, d_var, prob, map_comps, map_edges, costs_dict, inv_demand_dict, G, nodes = _construct_problem(edges, inv_edges)

    np.random.seed(seed)
    #draw random initial points
    r0_ = np.random.uniform(-10,10, (nsolutions, len(nodes)-1)) 

    r_tot = []
    #iteration over nsamples
    for i in range(nsolutions):
        r_k = []
        r_list = [s for s in r0_[i, :]]
        r_list.append(-np.sum(r_list))
        r.value = np.array(r_list)
        for j in range(max_iter):
            r_k.append(r.value)
            prob.solve(solver = cp.GUROBI)
            if prob.status!='optimal':
                print("iteration %d, status %s" %(i, prob.status))
            new_r = get_new_r(f_p, map_edges, nodes)
            r.value = new_r
        r_tot.append(r_k)

    return r_tot

def plot_results_run(r_tot, name = '', save = False):
    """
    Plot the results from running the algorithm

    Parameters
    ----------
    r_tot: list of lists
        each list in r_tot contains the evolution of r over iterations
    name: str
        name of the figure to be saved
    save: bool
        whether or not to save the different figures
    """
    n = len(r_tot[0][0])
    #Fig 1: plot the evolution of ri over iterations
    _, ax = plt.subplots(n, 1, figsize = (13, 10))
    for r_k in r_tot:
        r_arr = np.array(r_k)
        for i in range(n):
            ax[i].plot(r_arr[:,i], linewidth = 1)
    for i in range(n):
        ax[i].grid()
        ax[i].set_xlabel('$k$')
        ax[i].set_ylabel('Component %d' %i)
    if save:
        plt.savefig(os.path.join(IMAGES_PATH, name+"_r_k.png"), transparent = True)   

    #Fig 2: plot the evolution of ri over iterations, only for the last
    #nback iterations
    nback = 5
    _, ax = plt.subplots(n, 1, figsize = (13, 10))
    for r_k in r_tot:
        r_arr = np.array(r_k)
        for i in range(n):
            ax[i].plot(np.linspace(r_arr.shape[0]-nback, r_arr.shape[0],nback),r_arr[-nback:, i])
            # ax[i].plot(r_arr[-nback:,i], linewidth = 1)
    for i in range(n):
        ax[i].grid()
        ax[i].set_xlabel('$k$')
        ax[i].set_ylabel('Component %d' %i)
    if save:
        plt.savefig(os.path.join(IMAGES_PATH, name+"_r_k.png"), transparent = True)   

    #Fig 3: plot the evolution of the difference between r_k and r_k+1 
    #over iterations
    plt.figure()
    for r_k in r_tot:
        r_arr = np.array(r_k)
        diff = [np.linalg.norm(r_arr[i,:] - r_arr[i+1,:]) for i in range(r_arr.shape[0]-1)]
        plt.plot(diff)
    plt.grid()
    plt.xlabel('$k$')
    plt.ylabel('$\|r_k-r_{k+1}\|$')
    plt.yscale('log')
    if save: 
        plt.savefig(os.path.join(IMAGES_PATH, name+"_difference_rk.png"), transparent = True)   
    return


def sample_solutions(edges, inv_edges, nsamples=100, seed =0, name = '', save = False):
    """
    Draws a number nsamples of points and computes the mapping of those points via
    the optimization problem. 
    Useful diagnostic to see if map is contractive experimentally. 

    Parameters
    ----------
    name: str
        name under which to save figures if needed
    nsamples: int
        number of samples to draw
    seed: int 
        seed to initialize the random number generator
    """

    #build the optimization problem
    f_p, f_r, r, d_var, prob, map_comps, map_edges, costs_dict, inv_demand_dict, G, nodes = _construct_problem(edges, inv_edges)

    np.random.seed(seed)
    samples = np.random.uniform(-10,10, (nsamples, len(nodes)-1))

    r_ = []
    Tr_ = []
    #iteration over nsamples
    for i in range(nsamples):
        r_list = [s for s in samples[i,:]]
        r_list.append(-np.sum(r_list))
        r_.append(r_list)
        r.value = np.array(r_list)
        prob.solve(solver = cp.GUROBI)
        if prob.status!='optimal':
            print("iteration %d, status %s" %(i, prob.status))
        
        #computing the mapping
        Tr_.append(get_new_r(f_p, map_edges, nodes))

    dr= []
    dT = []
    for i in range(nsamples-1):
        dr.append(np.linalg.norm(
            np.array(r_[i]) - np.array(r_[i+1])
            ))
        dT.append(np.linalg.norm(
            np.array(Tr_[i]) - np.array(Tr_[i+1])
            ))

    #Fig 1: plot the distance of the mappings vs the original distances
    plt.figure(figsize=(10,6))
    plt.scatter(dr, np.divide(dT, dr), s = 10, marker = 'o')
    # plt.plot(dr, np.ones(len(dr)), 'r')
    plt.grid()
    plt.ylabel('$\|Tr_1 - Tr_2\|/\|r_1-r_2\|$')
    plt.xlabel('$\|r_1-r_2\|$')
    # plt.xlim([0, 20])
    # plt.ylim([0,1.3])
    if save:
        plt.savefig(os.path.join(IMAGES_PATH, name+".png"), transparent = True)    

    rat = np.divide(dT, dr)

    #give coordinates of points that have a ratio >1
    idx = np.where(rat>1)
    if idx[0].shape[0]>0:
        print("ratio larger than 1")
        for i in idx[0]: 
            print("----------------------")
            print("     values of ri: ", r_[i], r_[i+1])
            print("     values of Tr: ", Tr_[i], Tr_[i+1])
            print("     ratio: ", rat[i])
    return Tr_, r_, dT, dr

def get_new_r(f_p, map_edges, nodes):
    p_flows = f_p.value
    net_flows = dict()
    Tr = []
    for n in nodes:
        net_flows[n] = 0
        for e in map_edges:
            if e[0] == n:
                idx_list = map_edges[e]
                for ii in idx_list:
                    net_flows[n]+=p_flows[ii]
            elif e[1] ==n:
                idx_list = map_edges[e]
                for ii in idx_list:
                    net_flows[n]-=p_flows[ii]
        Tr.append(-net_flows[n])

    return Tr

def get_edge_flow(f_p, map_edges):
    edge_flow = dict()
    p_flows = f_p.value
    for e in map_edges:
        ii = map_edges[e]
        edge_flow[e] = 0
        for i in ii:
            edge_flow[e] += p_flows[i]
    return edge_flow