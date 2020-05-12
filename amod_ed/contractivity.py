import cvxpy as cp
from amod_ed.helpers_icu import BPR_int, BPR
import matplotlib.pyplot as plt
import numpy as np
import os

IMAGES_PATH = "/Users/lucasfuentes/ASL/Images"


def viz_costs(name, phi_p, phi_inv, k_p, k_inv, shift_inv, save = False):
    """
    Visualize the cost functions. 

    Parameters
    ----------
    name: str
        name of the figure if needs to be saved
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
    save: bool
        whether or not to save
    """
    edges = ["Edge $(1,2)$", "Edge $(2,1)$"]
    #plot paths for both directions
    x = np.linspace(0, 10, 100)
    _, ax = plt.subplots(2,1,figsize = (10, 10))

    #loop over both edges
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
    if save:
        plt.savefig(os.path.join(IMAGES_PATH, name+".png"), transparent = True) 

def viz_costs_natural(name, phi_p, phi_inv, k_p, k_inv, shift_inv, save = False):
    """
    We do the same as above except that the demand is linear, not BPR
    of the form 
    
        - phi_inv x + shift_inv

    Therefore, we ignore the values of k_inv. 

    """

    edges = ["Edge $(1,2)$", "Edge $(2,1)$"]
    #plot paths for both directions
    x = np.linspace(0, 10, 100)
    _, ax = plt.subplots(2,1,figsize = (10, 10))
    for i in range(2):
        cost = BPR(phi_p[i], x, k_p[i], beta = 4)
        inv_d = -phi_inv[i]*x+shift_inv[i]
        ax[i].plot(x, cost, label = 'Cost')
        ax[i].plot(x, inv_d, label = 'Inverse Demand')
        ax[i].grid()
        ax[i].set_xlabel("Total Flow")
        ax[i].set_ylabel("Cost")
        ax[i].legend()
        ax[i].set_title(edges[i])
        ax[i].set_ylim([cost[0]-10, inv_d[0]+10])
    if save:
        plt.savefig(os.path.join(IMAGES_PATH, name+".png"), transparent = True) 

def _construct_problem(phi_p, phi_inv, k_p, k_inv, shift_inv):
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

    #Define the optimization variables
    f_p = cp.Variable(2)
    f_r = cp.Variable(2)
    r = cp.Parameter()#in our setting, the rebalacing is a parameter

    constraints = [f_p>=0, f_r>=0, r == f_r[1]-f_r[0], f_p <= 10, f_r <=10]
    total_cost = 0
    #iterate over edges to build the total cost
    for i in range(2):
        #For each direction, there are two terms to the cost
        #1. the actual cost to the passenger
        #2. The inverse demand cost
        cost = BPR_int(phi_p[i], f_p[i]+f_r[i], k_p[i], beta = 4)
        inv_d = -BPR_int(phi_inv[i], f_p[i], k_inv[i], beta = 4) + shift_inv[i]*f_p[i]
        #The total cost is the sum of all those terms
        total_cost = total_cost + cost - inv_d
    #The objective is to minimize the total cost
    obj = cp.Minimize(total_cost)
    prob = cp.Problem(obj, constraints)
    return f_p, f_r, r, prob


def _construct_natural_problem(phi_p, phi_inv, k_p, k_inv, shift_inv):
    """
    Very similar to above, except that in this case the total cost 
    is only composed of the user cost (as we do not use the inverse demand 
    in the objective function but rather in the update). 
    """
    f_p = cp.Variable(2)
    f_r = cp.Variable(2)
    d = cp.Parameter(2, nonneg = True)
    
    constraints = [f_p>=0, f_r>=0, f_p + f_r == cp.max(d), f_p == d,]
    
    total_cost = 0
    for i in range(2):
        cost = BPR_int(phi_p[i], f_p[i]+f_r[i], k_p[i], beta = 4)
        total_cost = total_cost + cost
    obj = cp.Minimize(total_cost)
    prob = cp.Problem(obj, constraints)
    return f_p, f_r, d, prob


def _construct_initial_problem(phi_p, phi_inv, k_p, k_inv, shift_inv):
    """
    Build the optimization problem corresponding to the initial formulation. 
    In this case, we solve for everything at once. 

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

    #Define the optimization variables
    f_p = cp.Variable(2)
    f_r = cp.Variable(2)

    constraints = [f_p>=0, f_r>=0, f_p[0]-f_p[1] == f_r[1]-f_r[0], f_p <= 10, f_r <=10]
    total_cost = 0
    #iterate over edges to build the total cost
    for i in range(2):
        #For each direction, there are two terms to the cost
        #1. the actual cost to the passenger
        #2. The inverse demand cost
        cost = BPR_int(phi_p[i], f_p[i]+f_r[i], k_p[i], beta = 4)
        inv_d = -BPR_int(phi_inv[i], f_p[i], k_inv[i], beta = 4) + shift_inv[i]*f_p[i]
        #The total cost is the sum of all those terms
        total_cost = total_cost + cost - inv_d
    #The objective is to minimize the total cost
    obj = cp.Minimize(total_cost)
    prob = cp.Problem(obj, constraints)
    return f_p, f_r, prob

def run_algorithm(phi_p, phi_inv, k_p, k_inv, shift_inv, nsolutions = 5, seed =0, max_iter = 50):
    """
    Run the initial (not the natural) algorithm starting from a series of randomly
    sampled initial points

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
    f_p, _, r, prob = _construct_problem(phi_p, phi_inv, k_p, k_inv, shift_inv)
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
            r.value = f_p.value[0] - f_p.value[1]
        r_tot.append(r_k)

    return r_tot

def plot_results_run(r_tot, name, save = False):
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


def sample_solutions(name ,phi_p, phi_inv, k_p, k_inv, shift_inv, nsamples=100, seed =0, save = False):
    """
    Draws a number nsamples of points and computes the mapping of those points via
    the optimization problem. 
    Useful diagnostic to see if map is contractive experimentally. 

    Parameters
    ----------
    name: str
        name under which to save figures if needed
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
    nsamples: int
        number of samples to draw
    seed: int 
        seed to initialize the random number generator
    """

    #build the optimization problem
    f_p, _, r, prob = _construct_problem(phi_p, phi_inv, k_p, k_inv, shift_inv)

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
        Tr.append(f_p.value[0] - f_p.value[1])

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
    plt.ylim([0,1.3])
    if save:
        plt.savefig(os.path.join(IMAGES_PATH, name+".png"), transparent = True)    

    rat = np.divide(dT, dr)

    #give coordinates of points that have a ratio >1
    idx = np.where(rat>1)
    if idx[0].shape[0]>0:
        print("ratio larger than 1")
        for i in idx[0]: 
            print("     values of ri: ", samples[i], samples[i+1])
            print("     values of Tr: ", Tr[i], Tr[i+1])
    return dT, dr

def sample_natural_solutions(name, phi_p, phi_inv, k_p, k_inv, shift_inv, nsamples=100, seed =0, save =False):
    """
    Draws a number nsamples of points and computes the mapping of those points via
    the "natural" optimization problem. 
    Useful diagnostic to see if map is contractive experimentally. 

    Parameters
    ----------
    name: str
        name under which to save figures if needed
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
    nsamples: int
        number of samples to draw
    seed: int 
        seed to initialize the random number generator
    """

    #construct the natural problem
    f_p, f_r, d, prob = _construct_natural_problem(phi_p, phi_inv, k_p, k_inv, shift_inv)

    np.random.seed(seed)
    samples = np.random.uniform(0,10, (2, nsamples))
    Td = []

    #iterate over the number of samples
    for i in range(nsamples):
        d.value = samples[:,i]
        prob.solve(solver = cp.GUROBI)
        if prob.status!='optimal':
            print("iteration %d, status %s" %(i, prob.status))
        cost = BPR(phi_p, f_p.value+f_r.value, k_p, beta = 4)
        #apply the mapping of the demand
        #remember the demand is of the form 
        # - phi_inv x + shift_inv
        demand = - np.multiply(phi_inv, cost) + shift_inv
        Td.append(demand)

    dd= []
    dT = []
    #compute the different distances involved
    for i in range(nsamples-1):
        dd.append(np.linalg.norm(samples[:,i]-samples[:,i+1]))
        dT.append(np.linalg.norm(Td[i]-Td[i+1]))
    
    plt.figure(figsize=(10,6))
    plt.scatter(dd, np.divide(dT, dd), s = 10, marker = 'o')
    plt.grid()
    plt.ylabel('$\|Td_1 - Td_2\|/\|d_1-d_2\|$')
    plt.xlabel('$\|d_1-d_2\|$')
    if save:
        plt.savefig(os.path.join(IMAGES_PATH, name+".png"), transparent = True)    

    return dT, dd
