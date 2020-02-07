#This document contains the implementation of Kiril's algorithm to solve the assignment
#problem with rebalancers with graph extension, with Elastic Demand too

import numpy as np
import cvxpy as cp
import networkx as nx
from helpers_icu import Value_Total_Cost 
from routines_icu import update_costs, update_OD, update_capacities, AoN, estimate_ri_k, update_flows, init_flows, fixed_step
from result_analysis import check_flow_cons_at_OD_nodes

#TODO: think about the evolving bounds and the rebalancing smoothing
def solve(G_0, OD, edge_list, tol=10**-6, FW_tol=10**-6, max_iter=10**3):

    #Variables to store at each iterations
    i = 1
    G_ = []
    ri_ = []
    balance=[]
    #initialize certain values
    G_k = G_0
    ri_k, G_k = estimate_ri_k(G_k, ri_smoothing=False, a_k=0)
    n_iter_tot=0
    #Save the different variables
    G_.append(G_k)
    ri_.append(ri_k)

    compute = True
    try:
        while compute:
            print("##########################################")
            print("ITERATION #: ", i)
            # print("CURRENT RI_k")
            # print(ri_k)

            #TODO: maybe you do not have to go all the way in the computation
            #maybe you can introduce a decreasing tolerance over the different problems
            #like start with tol=1 and then divide it by two at every step
            #currently this is dealt with by a max iter number
            G_list, _, _, _ , n_iter= FW_graph_extension(
                G_k, OD, edge_list, ri_k, FW_tol,
                step='fixed', evolving_bounds=False, max_iter=max_iter)

            G_end = G_list[-1]

            #estimate #ri_k, update OD, assign, update costs
            ri_new, G_end = estimate_ri_k(G_end, ri_smoothing=False, a_k=0)

            #TODO: this is not a good measure I think
            if diff_ri(ri_k, ri_new) < tol:
                compute = False
                print("The rebalancing vector has reached a stationary point.")

            #update the values for the new iteration
            ri_k = ri_new
            # TODO: does it work if you actually keep the last version of G (as you solved it? )
            G_k = G_end
            balance.append(check_flow_cons_at_OD_nodes(G_k, OD))
            print("Balance norm at the end of iteration: ", np.linalg.norm(balance[-1]))
            #Save the different variables
            G_.append(G_k)
            ri_.append(ri_k)

            i += 1
            n_iter_tot+=n_iter
    except KeyboardInterrupt:
        print("Program interrupted by user -- Current data saved")
        return G_, ri_, i-1, n_iter_tot, np.array(balance)

    """
    Returns: 
    G_ : list of graphs
    ri_: list of dict()
    i: number of outer loop iterations, i.e. of ri update
    n_iter_tot: total number of iterations in the FW
    """
    return G_, ri_, i-1, n_iter_tot, np.array(balance)


def diff_ri(ri_k, ri_new):

    diff = []

    for n in ri_k.keys():
        diff.append(ri_k[n]-ri_new[n])
    diff = np.asarray(diff)
    return np.linalg.norm(diff)


def FW_graph_extension(G_0, OD, edge_list, ri_k, FW_tol=10**-6,
    step='line_search', evolving_bounds=True, max_iter=10**3):
    #ri_t are the estimate of ri at timestep k

    ###########################################
    # PARAMETERS
    ##########################################
    y_list = []
    opt_res = dict()
    opt_res['a_k'] = []
    opt_res['obj'] = []
    opt_res['dual_gap'] = []
    OD_list = []
    G_list = []
    G_list.append(G_0)
    G_k = G_0.copy()
    a_k = 1
    i = 1
    compute = True

    #################################
    # update the OD pairs and capacities
    #################################
    OD = update_OD(OD, ri_k, a_k, G_k, evolving_bounds)
    # print("CURRENT OD:", OD)
    #you update capacities because you have new values of ri_k
    G_k = update_capacities(G_k, ri_k)
    # we need to ensure that the information is passed on to the costs
    G_k = update_costs(G_k)

    ###################################
    # Reinitialize
    ###################################

    G_k = init_flows(G_k, OD)
    G_k = update_costs(G_k)

    ###################################
    # Solve for the given ri_k
    ###################################

    while compute:  

        #perform AON assignment
        y_k = AoN(G_k, OD)
        #TODO: enable line search
        if step == 'line_search':
            # a_k,obj_k=line_search(G_crt,y_k,edge_list)#include the fixed step size,
            print("not implemented")
            return
        elif step == 'fixed':
            a_k = fixed_step(i)
            obj_k=Value_Total_Cost(G_k)
        else:
            print("wrong optim step chosen")
            return

        #compute the duality gap
        duality_gap = compute_duality_gap(G_k, y_k)
        if duality_gap < FW_tol or i >= max_iter:
            #here we put a limit on the number of computations as the problem is likely to change
            #however we do not do it in the main loop!
            #TODO: I believe this stopping criterion scheme is really not reliable, as it does not
            #guarantee anything...
            compute = False
            if duality_gap < FW_tol:
                print("     FW solved to tol")
            else:
                print("     Max inner iterations reached")
            print("     Number of inner loop iterations: ", i)

        #update the flows
        G_k = update_flows(G_k, y_k, a_k, edge_list)
        G_k = update_costs(G_k)

        #save for analyses
        opt_res['obj'].append(obj_k)
        opt_res['a_k'].append(a_k)
        opt_res['dual_gap'].append(duality_gap)
        G_list.append(G_k)
        y_list.append(y_k)
        OD_list.append(OD.copy())
        i += 1
    return G_list, y_list, opt_res, OD_list, i-1


def compute_duality_gap(G_k, y_k):
    #G_k : version of the graph (and therefore the flows) at iteration k
    #y_k : minimizer of linearized problem at iteration k
    d_gap = 0
    for e in G_k.edges():
        for flag in ['f_m', 'f_r']:
            x_k_ij = G_k[e[0]][e[1]][flag]
            c_k_ij = G_k[e[0]][e[1]]['cost']
            y_k_ij = y_k[(e[0], e[1]), flag]

            # print("EDGE: ", e, " | flag: " , flag)
            # print("x: ", x_k_ij, " | y: ", y_k_ij, " | c: ", c_k_ij)
            # print("update: ", (x_k_ij-y_k_ij)*c_k_ij)
            d_gap += (x_k_ij-y_k_ij)*c_k_ij

    return d_gap








###########################################################################
# Original version of functions in FW_icu
#
# Those functions are already present in the FW_icu file. We need another
# version because here we want to update all flows at the same time
###########################################################################

#TODO: consolidate them all in a single file, with flexible versions


#TODO: Deprecate

#main problem in the below: the flows are updated twice in the "same" direction
#as the updates occur in place I think
#plus, there is absolutely no need to compute the flows here... 

# def fixed_step(G, y_k, edge_list, k):
#     """
#     """

#     #here we actually update both the rebalancing and the passenger flows
#     #we need to take both into account because we are solving for a fixed value of the ri's

#     # the update step is very important, if it is too large, then we lose the progress we made
#     gamma = 2/(k+2)
#     #currently I put k**1.5 because I want to accelerate the computation slightly

#     for i in range(len(edge_list)):
#         e = edge_list[i]
#         for flag in ['f_m', 'f_r']:
#             x_k_e = G[e[0]][e[1]][flag]  # retrieve the flow
#             # retrieve the flow from the manual assignment
#             y_k_e = y_k[(e[0], e[1]), flag]
#             G[e[0]][e[1]][flag] = (1-gamma)*x_k_e+gamma*y_k_e

#     return gamma, Value_Total_Cost(G)


# def update_flows(G, y_k, a_k, edge_list):

#     for i in range(len(edge_list)):
#         e = edge_list[i]
#         for flag in ['f_m', 'f_r']:
#             x_k_e = G[e[0]][e[1]][flag]  # retrieve the flow
#             # retrieve the flow from the manual assignment
#             y_k_e = y_k[(e[0], e[1]), flag]
#             G[e[0]][e[1]][flag] = (1-a_k)*x_k_e + a_k * y_k_e
#     return G


# def init_flows(G, OD):
#     #TODO: deal with conflicting version!!
#     #Initiliaze the flows, with a feasible solution
#     #still unclear whether we need to initialize the rebalancers appropriately too

#     #reinitialize the flows to zero
#     for e in G.edges():
#         for flag in ['f_r', 'f_m']:
#             G[e[0]][e[1]][flag] = 0

#     for (o, d) in OD.keys():
#         path = nx.shortest_path(G, source=o, target=d, weight='cost')
#         for i in range(len(path)-1):
#             if d == 'R':
#                 G[path[i]][path[i+1]]['f_r'] += OD[o, d]
#             else:
#                 G[path[i]][path[i+1]]['f_m'] += OD[o, d]
#     return G
