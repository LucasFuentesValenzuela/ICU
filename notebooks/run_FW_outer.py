import os
from graph import construct_graph, get_edge_list
from FW_OuterUpdate import solve
import pickle
import argparse
from result_analysis import plot_ri, print_balance

def main():
    """
    Runs the algorithm on the datafiles given as arguments
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", help="path to data directory.",
                        type=str)
            
    parser.add_argument("-L", help="Length of rebalancing edge.",
                        type=int, default=10**4)
    
    parser.add_argument("-ni", help="Maximum number of iterations in the inner FW loop.",
                        type=int, default=5000)

    parser.add_argument("-no", help="Maximum number of iterations in the outer loop.",
                        type=int, default=50) 
    parser.add_argument("-sd", help="Saving Directory for this experiment", 
                        type=str, default='.')

    parser.add_argument("-ev", help="Evolving bounds", type=int, default=0)

    parser.add_argument("-sc", help="Stopping criterion", type=str, default="rp")

    parser.add_argument("-fu", help='Update of the step to avoid overshooting', type=int, default=1)

    args = parser.parse_args()
    path = os.path.join('Data/',args.p)
    save_dir=os.path.join(path,'outputs',args.sd)
    L_r=args.L
    ni=args.ni
    no=args.no
    ev=args.ev
    sc=args.sc
    fu=args.fu

    if ev==1:
        evolving_bounds=True
    else:
        evolving_bounds=False
    if fu==1:
        fixed_update=True
    else:
        fixed_update=False

    if not os.path.exists(path):
        print("Path does not exist")
        return
        
    if sc=='rp':
        stopping_criterion="relative_progress"
    elif sc=='dg':
        stopping_criterion='duality_gap'
    else:
        print("wrong stopping criterion chosen")
        return


    print("---------------------------------------------------")
    print("Solving with the following parameters:")
    print("     path: ", path)
    print("     L_r: ", L_r)
    print("     ni: ", ni)
    print("     no: ", no)
    print("     ev: ", evolving_bounds)
    print("     sc: ", stopping_criterion)
    print("     fu: ", fixed_update)

    G_0, OD = construct_graph(path, L_rebalancing_edge = L_r)
    edge_list = get_edge_list(G_0)
    G_FW, ri_FW, n_outer, n_inner, balance, opt_res, OD_list, balance_list = solve(
        G_0.copy(), OD.copy(), edge_list, tol=10**-3, FW_tol=0, 
        max_iter_outer= no, max_iter= ni, evolving_bounds=evolving_bounds,
        stopping_criterion = stopping_criterion, update=fixed_update,
        ri_smoothing=False)
    
    #Save the different documents
    #TODO: save the parameters too!! 
    to_save=[G_FW, OD, ri_FW, n_outer, n_inner, balance, opt_res, OD_list, balance_list]

    print('Saving to external .pkl file')
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    with open(os.path.join(save_dir,'output'+'_L_'+str(L_r)+'_ni_'+str(ni)+'_no_'+str(no)+'_ev_'+str(ev)+
    '_'+stopping_criterion+'_fu_'+str(fixed_update)+'.pkl'), 'wb') as f:
            pickle.dump(to_save, f)

    return


if __name__ == "__main__":
    main()
