# This document contains the implementation of Kiril's algorithm to solve the assignment
# problem with rebalancers with graph extension, with Elastic Demand too

import numpy as np
import networkx as nx
from amod_ed.helpers_icu import Value_Total_Cost
from amod_ed.routines_icu import update_costs, update_OD, update_capacities, AoN, estimate_ri_k, update_flows, fixed_step, line_search
from amod_ed.routines_icu import backtracking_line_search
from amod_ed.routines_icu import check_flow_cons_at_OD_nodes
from amod_ed.flows_init import initialize_flows

# TODO: think about the evolving bounds and the rebalancing smoothing
# TODO: update stopping criterion (should be based on balance norm I think)
# TODO: check duality gap
# TODO: make line search work
# TODO: make what we save consistent (same length, etc.)
# TODO: decide what to do with evolving tol


def solve(
        G_0, OD, edge_list, tol=10**-6, FW_tol=10**-6, max_iter_outer=50,
        max_iter=10**3, evolving_bounds=True,
        stopping_criterion='relative_progress', update_factor=1.1,
        ri_smoothing=True):

    # Variables to store at each iterations
    i = 1
    G_ = []
    ri_ = []
    balance = []
    balance_list = []  # list of list of balances
    opt_res_ = []
    OD_list = []

    # initialize certain values
    G_k = G_0
    ri_k, G_k = estimate_ri_k(G_k, ri_smoothing=False, a_k=0, ri_prev=[])
    balance_k = check_flow_cons_at_OD_nodes(G_k, OD)
    balance_k = np.ones(balance_k.shape)
    n_iter_tot = []

    # Save the different variables
    G_.append(G_k)
    ri_.append(ri_k)
    OD_list.append(OD)
    balance.append(balance_k)
    # opt_res_.append([])

    FW_tol_k = FW_tol
    compute = True

    # True when we reached a stage where we want to keep the previous value of the step
    # useful only if we use fixed step
    # TODO: implement that when using fixed step.
    continuous_step = False #I think you can enable it relatively early on
    i_offset = 0

    try:
        while compute:

            # TODO: make this the standard
            if i > 10:
                continuous_step = True

            print("##########################################")
            print("ITERATION #: ", i)
            print("Current FW tol: ", FW_tol_k)
            print("Current max ni: ", max_iter)
            print("i_offset : ", i_offset)

            G_list, _, opt_res_k, OD_list_k, n_iter, balance_ = FW_graph_extension(
                G_k.copy(), OD.copy(), edge_list, ri_k, FW_tol=FW_tol_k,
                step='line_search', evolving_bounds=evolving_bounds, max_iter=max_iter,
                stopping_criterion=stopping_criterion, update_factor=update_factor, i_offset=i_offset)

            # this is a good choice only if you have monotonous decrease
            G_end = G_list[-1].copy()

            a_k_smoothing = 0 
            # if i > 10:
                # ri_smoothing = True
            print("ri_smoothing: ", ri_smoothing)

            ri_new, G_end = estimate_ri_k(
                G_end.copy(), ri_smoothing=ri_smoothing, a_k=a_k_smoothing, ri_prev=ri_k)

            balance_new = check_flow_cons_at_OD_nodes(G_end.copy(), OD)
            balance_norm = np.linalg.norm(balance_new)

            # previous criterion based on the difference of balanced
            # diff_balance = np.linalg.norm(balance_new-balance_k)

            if balance_norm < tol:
                compute = False
                print("The balance vector has reached the tol.")
            elif i >= max_iter_outer:
                compute = False
                print("Maximum number of outer iterations reached")

            ri_k = ri_new

            G_k = G_end.copy()
            balance_k = balance_new

            print("Balance norm at the end of iteration: ",
                  np.around(balance_norm, 2))

            # Save the different variables
            G_.append(G_list)
            balance.append(balance_k)
            opt_res_.append(opt_res_k)
            n_iter_tot.append(n_iter)
            OD_list.append(OD_list_k)
            balance_list.append(balance_)
            ri_.append(ri_k)

            # TODO: if you want to do evolving tol
            # FW_tol_k = np.maximum(FW_tol, FW_tol_k/2)

            # TODO: implement this for real
            if continuous_step:
                i_offset += n_iter
                # i_offset = i_offset ** 1.5
            i += 1
    except KeyboardInterrupt:
        print("Program interrupted by user -- Current data saved")
        return G_, ri_, i-1, n_iter_tot, np.array(balance), opt_res_, OD_list, balance_list

    return G_, ri_, i-1, n_iter_tot, np.array(balance), opt_res_, OD_list, balance_list

# TODO: consolidate what we think about diff ri, diff f_r, balance... as stopping crit. 
# def diff_ri(ri_k, ri_new):

#     diff = []

#     for n in ri_k.keys():
#         diff.append(ri_k[n]-ri_new[n])
#     diff = np.asarray(diff)
#     return np.linalg.norm(diff)


def FW_graph_extension(G_0, OD, edge_list, ri_k, FW_tol=10**-6,
                       step='fixed', evolving_bounds=False, max_iter=10**3,
                       stopping_criterion='relative_progress',
                       update_factor=1.1, i_offset=0):

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
    balance_list = []

    G_k = G_0.copy()

    a_k = 1
    i = 1

    compute = True

    #################################
    # update the OD pairs and capacities
    #################################
    OD = update_OD(OD, ri_k, G_k, evolving_bounds)
    # you update capacities because you have new values of ri_k
    G_k = update_capacities(G_k, ri_k)
    # we need to ensure that the information is passed on to the costs
    G_k = update_costs(G_k)

    ###################################
    # Reinitialize
    ###################################

    # replace None by G_0 if want proper init
    G_k = initialize_flows(G_k, G_0, ri_k, OD) #we do not use init with NN just to use line search and see

    G_k = update_costs(G_k)

    obj_k, G_k = Value_Total_Cost(G_k)

    opt_res['obj'].append(obj_k)
    G_list.append(G_k.copy())
    ###################################
    # Solve for the given ri_k
    ###################################

    while compute:
        i_step = i + i_offset

        # compute the balance
        balance_new = check_flow_cons_at_OD_nodes(G_k.copy(), OD)
        balance_norm = np.linalg.norm(
            balance_new)/np.sqrt(balance_new.shape[0])
        balance_list.append(balance_norm)

        ############################################
        # ASSIGNMENT
        ############################################
        y_k = AoN(G_k.copy(), OD)

        if step == 'line_search':
            # TODO: make this work
            # a_k = line_search(G_k, y_k, edge_list)
            a_k = backtracking_line_search(G_k, y_k, a = 0.1, b = .8)

        elif step == 'fixed':
            # TODO: decide whether to scrap the update factor
            if isinstance(update_factor, float):
                a_k = fixed_step(i, y_k, G_k, update=True,
                                 update_factor=update_factor)
            else:
                a_k = fixed_step(i_step, y_k, G_k, update=False)
        else:
            print("wrong optim step chosen")
            return

        ######################################
        # Stopping criterion
        #######################################
        opt_res, compute = compute_stopping_criterion(
            stopping_criterion, i, G_k, y_k, opt_res,
            FW_tol, compute, max_iter, n_rolling=5)
        
        if not compute:
            continue

        # update the flows
        G_k = update_flows(G_k, y_k, a_k, edge_list)
        G_k = update_costs(G_k)

        if evolving_bounds:
            # you need to update the OD (evolving bounds)
            OD = update_OD(OD, ri_k, G_k, evolving_bounds)

        # save for analyses
        obj_k, G_k = Value_Total_Cost(G_k.copy())
        opt_res['obj'].append(obj_k)
        opt_res['a_k'].append(a_k)
        G_list.append(G_k.copy())
        y_list.append(y_k.copy())
        OD_list.append(OD.copy())
        i += 1

    return G_list, y_list, opt_res, OD_list, i, balance_list


def compute_duality_gap(G_k, y_k):
    # G_k : version of the graph (and therefore the flows) at iteration k
    # y_k : minimizer of linearized problem at iteration k
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


def compute_stopping_criterion(
        stopping_criterion, i, G_k, y_k, opt_res,
        FW_tol, compute, max_iter, n_rolling=5):

    if stopping_criterion == 'duality_gap' and i > 1:
        # compute the duality gap
        duality_gap = compute_duality_gap(G_k, y_k)
        opt_res['stop'].append(duality_gap)
        if duality_gap < FW_tol or i >= max_iter:
            # here we put a limit on the number of computations as the problem is likely to change
            # however we do not do it in the main loop!
            compute = False
            if duality_gap < FW_tol:
                print("     FW solved to tol")
            else:
                print("     Max inner iterations reached")
            print("     Number of inner loop iterations: ", i)

    elif stopping_criterion == 'relative_progress' and i > 1:
        obj_prev = opt_res['obj'][-2]
        obj_k = opt_res['obj'][-1]
        rel_progress = abs(obj_prev-obj_k)/obj_prev
        opt_res['stop'].append(rel_progress)

        #We evaluate it over a rolling mean to make sure we just did not hit a random lucky spot
        if len(opt_res['stop']) > n_rolling and np.mean(opt_res['stop'][-n_rolling:]) < FW_tol:
            compute = False
            print("     Number of inner loop iterations: ", i-1)
            print("     FW solved to tol")

        elif i > max_iter:
            compute = False
            print("    Max inner iterations reached")
            print("     Number of inner loop iterations: ", i-1)

    return opt_res, compute
