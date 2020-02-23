#This document contains the implementation of Kiril's algorithm to solve the assignment
#problem with rebalancers with graph extension, with Elastic Demand too

import numpy as np
# import cvxpy as cp
import networkx as nx
from helpers_icu import Value_Total_Cost 
from routines_icu import update_costs, update_OD, update_capacities, AoN, estimate_ri_k, update_flows, init_flows, fixed_step
from result_analysis import check_flow_cons_at_OD_nodes

#TODO: think about the evolving bounds and the rebalancing smoothing
#TODO: update stopping criterion (should be based on balance norm I think)
#TODO: check duality gap

def solve(
    G_0, OD, edge_list, tol=10**-6, FW_tol=10**-6, max_iter_outer = 50, 
    max_iter=10**3, evolving_bounds=True, 
    stopping_criterion = 'relative_progress', update=True,
    ri_smoothing=True):

    #Variables to store at each iterations
    i = 1
    G_ = []
    ri_ = []
    balance=[]
    balance_list=[]#list of list of balances
    opt_res_=[]
    OD_list=[]
    #initialize certain values
    G_k = G_0
    ri_k, G_k = estimate_ri_k(G_k, ri_smoothing=False, a_k=0, ri_prev=[])
    balance_k=check_flow_cons_at_OD_nodes(G_k, OD)
    balance_k=np.ones(balance_k.shape)
    n_iter_tot=[]
    #Save the different variables
    G_.append(G_k)
    ri_.append(ri_k)
    OD_list.append(OD)
    balance.append(balance_k)
    # opt_res_.append([])

    # FW_tol_k=.1
    FW_tol_k=FW_tol

    compute = True
    try:
        while compute:
            print("##########################################")
            print("ITERATION #: ", i)
            print("Current FW tol: ", FW_tol_k)
            # print(ri_k)

            #TODO: maybe you do not have to go all the way in the computation
            #maybe you can introduce a decreasing tolerance over the different problems
            #like start with tol=1 and then divide it by two at every step
            #currently this is dealt with by a max iter number
            G_list, _, opt_res_k, OD_list_k , n_iter, balance_ = FW_graph_extension(
                G_k.copy(), OD.copy(), edge_list, ri_k, FW_tol_k,
                step='fixed', evolving_bounds=evolving_bounds, max_iter=max_iter,
                stopping_criterion = stopping_criterion, update=update)

            G_end = G_list[-1]#this is a good choice only if you have monotonous decrease

            #estimate #ri_k, update OD, assign, update costs
            ri_new, G_end = estimate_ri_k(G_end, ri_smoothing=ri_smoothing, a_k=1/2, ri_prev=ri_k)

            # for n in G_end.nodes():
            #     print(n, " ri: ", ri_new[n])
            #     print(n, " G_ri: ", G_end.nodes[n]["ri"])

            balance_new=check_flow_cons_at_OD_nodes(G_end, OD)
            balance_norm=np.linalg.norm(balance_new)
            diff_balance=np.linalg.norm(balance_new-balance_k)

            #New stopping criterion
            #Now based on the balance
            #We might still get stuck on an local optimum? 
            #Hence, maybe include a condition on the actual norm of the balance vector
            #for instance np.linalg.norm(balance_k)<eps

            #TODO: make sure we do not have to check that across ALL nodes
            # if diff_ri(ri_k, ri_new) < tol or i>=max_iter_outer:
            #     compute=False

            if diff_balance<tol:
                compute=False
                print("The balance vector has reached a stationary point")
            elif i>=max_iter_outer:
                compute = False
                print("Maximum number of outer iterations reached")

            #update the values for the new iteration
            ri_k = ri_new
            #some kind of rebalancer smoothing
            #TODO: integrate in the estimate ri_k routine
            r=dict()
            for n in ri_k.keys():
                r[n]=[]
            for n in r.keys():
                for ri in ri_:
                    r[n].append(ri[n])
                    ind=np.minimum(4,len(r[n]))
                    avg_n=np.mean(r[n][-ind:])
                    beta=2/(i+2)    
                    ri_k[n]=(beta)*ri_k[n]+(1-beta)*avg_n
            ri_.append(ri_k)

            # TODO: does it work if you actually keep the last version of G (as you solved it? )
            G_k = G_end
            balance_k=balance_new

            #the below might be completely wrong
            #indeed, OD is not complete here as it does not contain the rebalancing 
            #therefore, surely it is smaller than the actual OD in the dataset
            print("Balance norm at the end of iteration: ", np.around(balance_norm,2))

            #Save the different variables
            G_.append(G_list)
            
            balance.append(balance_k)
            opt_res_.append(opt_res_k)
            n_iter_tot.append(n_iter)
            OD_list.append(OD_list_k)
            balance_list.append(balance_)

            FW_tol_k=np.maximum(FW_tol,FW_tol_k/2)

            i += 1
    except KeyboardInterrupt:
        print("Program interrupted by user -- Current data saved")
        return G_, ri_, i-1, n_iter_tot, np.array(balance), opt_res_, OD_list, balance_list

    """
    Returns: 
    G_ : list of graphs
    ri_: list of dict()
    i: number of outer loop iterations, i.e. of ri update
    n_iter_tot: total number of iterations in the FW
    """
    return G_, ri_, i-1, n_iter_tot, np.array(balance), opt_res_, OD_list, balance_list


def diff_ri(ri_k, ri_new):

    diff = []

    for n in ri_k.keys():
        diff.append(ri_k[n]-ri_new[n])
    diff = np.asarray(diff)
    return np.linalg.norm(diff)


def FW_graph_extension(G_0, OD, edge_list, ri_k, FW_tol=10**-6,
    step='fixed', evolving_bounds=False, max_iter=10**3, 
    stopping_criterion = 'relative_progress', update=True):
    #ri_t are the estimate of ri at timestep k

    ###########################################
    # PARAMETERS
    ##########################################
    y_list = []
    opt_res = dict()
    opt_res['a_k'] = []
    opt_res['obj'] = []
    opt_res['stop'] = []
    OD_list = []
    G_list = []
    balance_list=[]

    G_k = G_0.copy()

    a_k = 1
    i = 1

    compute = True

    #################################
    # update the OD pairs and capacities
    #################################
    OD = update_OD(OD, ri_k, G_k, evolving_bounds)
    #you update capacities because you have new values of ri_k
    G_k = update_capacities(G_k, ri_k)
    # we need to ensure that the information is passed on to the costs
    G_k = update_costs(G_k)

    ###################################
    # Reinitialize
    ###################################

    G_k = init_flows(G_k, OD) #I think this is stupid: you lose so much. You should make the current state feasible (i.e. take away some flow but keep it feasible)
    G_k = update_costs(G_k)

    G_list.append(G_k.copy()) 
    obj_k, G_k = Value_Total_Cost(G_k)
    opt_res['obj'].append(obj_k)
    ###################################
    # Solve for the given ri_k
    ###################################

    while compute:  
        #TODO: figure out what the best order between assignment and stopping criterion is
        ############################################
        # ASSIGNMENT
        ############################################
        #perform AON assignment
        y_k = AoN(G_k, OD)

        #TODO: enable line search
        if step == 'line_search':
            # a_k,obj_k=line_search(G_crt,y_k,edge_list)#include the fixed step size,
            print("not implemented")
            return
        elif step == 'fixed':
            a_k = fixed_step(i, y_k, G_k, update=update)
        else:
            print("wrong optim step chosen")
            return
        ############################################
        ######################################
        # Stopping criterion
        #######################################
        n_rolling=5

        if stopping_criterion == 'duality_gap' and i>1:
            #compute the duality gap
            duality_gap = compute_duality_gap(G_k, y_k)
            opt_res['stop'].append(duality_gap)
            if duality_gap < FW_tol or i >= max_iter:
                #here we put a limit on the number of computations as the problem is likely to change
                #however we do not do it in the main loop!
                compute = False
                if duality_gap < FW_tol:
                    print("     FW solved to tol")
                else:
                    print("     Max inner iterations reached")
                print("     Number of inner loop iterations: ", i)
                continue
        elif stopping_criterion == 'relative_progress' and i>1:
            obj_prev = opt_res['obj'][-2]
            obj_k = opt_res['obj'][-1]
            rel_progress= abs(obj_prev-obj_k)/obj_prev
            opt_res['stop'].append(rel_progress)
            if len(opt_res['stop'])>n_rolling and np.mean(opt_res['stop'][-n_rolling:]) < FW_tol:
                compute = False
                if rel_progress < FW_tol:
                    print("     Number of inner loop iterations: ", i)
                    print("     FW solved to tol")
                    continue
            elif i >= max_iter:
                compute=False
                print("    Max inner iterations reached")
                print("     Number of inner loop iterations: ", i)
                continue
        
        
        #update the flows
        G_k = update_flows(G_k, y_k, a_k, edge_list)
        G_k = update_costs(G_k)

        if evolving_bounds:
            OD = update_OD(OD, ri_k, G_k, evolving_bounds) #you need to update the OD (evolving bounds)


        balance_new=check_flow_cons_at_OD_nodes(G_k, OD)
        balance_norm=np.linalg.norm(balance_new)/np.sqrt(balance_new.shape[0])
        balance_list.append(balance_norm)
        #save for analyses
        obj_k, G_k = Value_Total_Cost(G_k.copy())
        opt_res['obj'].append(obj_k)
        opt_res['a_k'].append(a_k)
        G_list.append(G_k.copy())
        y_list.append(y_k.copy())
        OD_list.append(OD.copy())
        i += 1
    return G_list, y_list, opt_res, OD_list, i-1, balance_list


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