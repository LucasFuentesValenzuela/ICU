import cvxpy as cp
from amod_ed.helpers_icu import BPR_int, BPR
import matplotlib.pyplot as plt
import numpy as np
import os


"""
A more complex and complete version of the contractivity framework
"""
IMAGES_PATH = "/Users/lucasfuentes/ASL/Images"


def viz_costs(edges, inv_edges, name = '', save = False):
    """
    Visualize the cost functions. 

    Parameters
    ----------
    save: bool
        whether or not to save
    """

    #plot paths for both directions
    x = np.linspace(0, 10, 100)
    _, ax = plt.subplots(len(edges),1,figsize = (10, 5*len(edges)))

    #loop over both edges
    for i in range(len(edges)):
        o = edges.loc[i, 'origin']
        d = edges.loc[i, 'destination']
        
        phi_p = edges.loc[i, 'phi']
        k_p = edges.loc[i, 'k']

        ind_inv = inv_edges.index[np.logical_and(
            inv_edges['origin'] == o, inv_edges['destination'] == d)][0]
        phi_inv = inv_edges.loc[ind_inv, 'phi']
        k_inv = inv_edges.loc[ind_inv, 'k']
        shift_inv = inv_edges.loc[ind_inv, 'shift']

        cost = BPR(phi_p, x, k_p, beta = 4)
        inv_d = -BPR(phi_inv, x, k_inv, beta = 4) + shift_inv

        ax[i].plot(x, cost, label = 'Cost')
        ax[i].plot(x, inv_d, label = 'Inverse Demand')
        ax[i].grid()
        ax[i].set_xlabel("Total Flow")
        ax[i].set_ylabel("Cost")
        ax[i].legend()
        ax[i].set_title('Edge (%d, %d)'%(o,d))
        ax[i].set_ylim([cost[0]-10, inv_d[0]+10])
    if save:
        plt.savefig(os.path.join(IMAGES_PATH, name+".png"), transparent = True) 

def _construct_problem(edges, inv_edges):
    """
    Build the optimization problem based on the shape of the cost functions. 

    Parameters
    ----------
    phi_p: list
        list of floats, containing the value of phi for the passengers for each edge
    phi_inv: list
        list of floats, containing the value of phi for the inverse demand functions
    k_p: list
        list of floats, containing the value of kappa for each edge
    k_inv: list
        list of floats, containing the value of kappa for the inverse demand edges
    shift_inv: list
        list of floats, containing the value of the inverse demand shifts 
    """
    #get the edge list
    edges_ = [(edges.loc[i, 'origin'], edges.loc[i, 'destination']) for i in edges.index]
    edge_list = []
    for e in edges_:
        if e not in edge_list:
            edge_list.append(e)

    #identifies each edge with a "direction number"
    n_edges =np.zeros(len(edges_))

    for j in range(len(edges)):
        if edges_[j] == edge_list[1]:
            n_edges[j]=1


    #Define the optimization variables
    f_p = cp.Variable(len(edges))
    f_r = cp.Variable(len(edges))
    r = cp.Parameter()#in our setting, the rebalacing is a parameter

    constraints = [f_p>=0, f_r>=0, r == (2*n_edges -1)*f_r, f_p <= 10, f_r <=10]
    total_cost = 0
    #iterate over edges to build the total cost
    #build the edge costs
    for i in range(len(edges)):
        #For each direction, there are two terms to the cost
        #1. the actual cost to the passenger
        #2. The inverse demand cost
        phi_p = edges.loc[i, 'phi']
        k_p = edges.loc[i,'k']
        cost = BPR_int(phi_p, f_p[i]+f_r[i], k_p, beta = 4)
        #The total cost is the sum of all those terms
        total_cost = total_cost + cost
    for i in range(len(inv_edges)):

        flow_p = np.sum((n_edges == i) * f_p)

        k_inv = inv_edges.loc[i, 'k']
        phi_inv = inv_edges.loc[i, 'phi']
        shift = inv_edges.loc[i, 'shift']
        inv_d = -BPR_int(phi_inv, flow_p, k_inv, beta = 4) + shift*flow_p
        total_cost = total_cost - inv_d

    #The objective is to minimize the total cost
    obj = cp.Minimize(total_cost)
    prob = cp.Problem(obj, constraints)
    return f_p, f_r, r, prob, n_edges





def run_algorithm(edges, inv_edges, nsolutions = 5, seed =0, max_iter = 50):
    """
    Run the initial (not the natural) algorithm starting from a series of randomly
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

    #construct the optimization problem
    f_p, _, r, prob, n_edges = _construct_problem(edges, inv_edges)
    #initialize the seed
    np.random.seed(seed)

    #draw random initial points
    r0_ = np.random.uniform(-10,10, nsolutions) 
    r_tot =[]
    #iterate for each initial point
    for i in range(nsolutions):
        r_k=[]
        r.value = r0_[i]
        for j in range(max_iter):
            r_k.append(r.value)
            prob.solve(solver=cp.GUROBI)
            net_flow = [np.sum((n_edges == i)*f_p.value) for i in np.unique(n_edges)]
            # balance = [np.sum((n_edges == i)*f_p.value) for i in np.unique(n_edges)]
            r.value = net_flow[0]-net_flow[1]
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

    #Fig 1: plot the evolution of ri over iterations
    plt.figure()
    for r_k in r_tot:
        plt.plot(r_k)
    plt.grid()
    plt.xlabel('$k$')
    plt.ylabel('$r_k$')
    if save:
        plt.savefig(os.path.join(IMAGES_PATH, name+"_r_k.png"), transparent = True)   

    #Fig 2: plot the evolution of ri over iterations, only for the last
    #nback iterations
    nback = 10
    plt.figure()
    for r_k in r_tot:
        plt.plot(np.linspace(len(r_k)-nback, len(r_k),nback),r_k[-nback:])
    plt.grid()
    plt.xlabel('$k$')
    plt.ylabel('$r_k$')
    if save:
        plt.savefig(os.path.join(IMAGES_PATH, name+"_r_k_zoom.png"), transparent = True) 

    #Fig 3: plot the evolution of the difference between r_k and r_k+1 
    #over iterations
    plt.figure()
    for r_k in r_tot:
        diff = [np.abs(r_k[i] - r_k[i+1]) for i in range(len(r_k)-1)]
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
    f_p, f_r, r, prob, n_edges = _construct_problem(edges, inv_edges)

    np.random.seed(seed)
    samples = np.random.uniform(-10,10, nsamples)

    #Tr denotes the mapping of T acting on r
    Tr = []
    #iteration over nsamples
    for i in range(nsamples):
        r.value = samples[i]
        prob.solve(solver = cp.GUROBI)
        if prob.status!='optimal':
            print("iteration %d, status %s" %(i, prob.status))
        net_flow = [np.sum((n_edges == i)*f_p.value) for i in np.unique(n_edges)]
        Tr.append(net_flow[0]-net_flow[1])

    dr= []
    dT = []
    for i in range(nsamples-1):
        dr.append(np.abs(samples[i]-samples[i+1]))
        dT.append(np.abs(Tr[i]-Tr[i+1]))

    #Fig 1: plot the distance of the mappings vs the original distances
    plt.figure(figsize=(10,6))
    plt.scatter(dr, np.divide(dT, dr), s = 10, marker = 'o')
    # plt.plot(dr, np.ones(len(dr)), 'r')
    plt.grid()
    plt.ylabel('$\|Tr_1 - Tr_2\|/\|r_1-r_2\|$')
    plt.xlabel('$\|r_1-r_2\|$')
    plt.xlim([0, 20])
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
            print("     values of ri: ", samples[i], samples[i+1])
            print("     values of Tr: ", Tr[i], Tr[i+1])
            print("     ratio: ", rat[i])
    return dT, dr
