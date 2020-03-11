import os
import pickle
import argparse
from amod_ed.result_analysis import plot_ri, print_balance
from amod_ed.FW_OuterUpdate import solve
from amod_ed.graph import construct_graph, get_edge_list


DATA_PATH = '/Users/lucasfuentes/ASL/ICU/notebooks/Data'


#TODO: minimize the amount of data saved when running experiments
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

    parser.add_argument("-fu", help='Update of the step to avoid overshooting', type=float, default=1.1)

    args = parser.parse_args()
    path = os.path.join(DATA_PATH,args.p)
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
    if fu==0:
        fixed_update=False
        update_factor=False
    else:
        fixed_update=True
        update_factor=float(fu)

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


    ri_smoothing=False
    FW_tol=10**-6
    tol=0

    G_0, OD = construct_graph(path, L_rebalancing_edge = L_r)
    edge_list = get_edge_list(G_0)
    G_FW, ri_FW, n_outer, n_inner, balance, opt_res, OD_list, balance_list = solve(
        G_0.copy(), OD.copy(), edge_list, tol=tol, FW_tol=FW_tol, 
        max_iter_outer= no, max_iter= ni, evolving_bounds=evolving_bounds,
        stopping_criterion = stopping_criterion, update_factor=update_factor,
        ri_smoothing=ri_smoothing)


    params=dict()
    params['L']=L_r
    params['ni']=ni
    params['no']=no
    params['ev']=evolving_bounds
    params['sc']=stopping_criterion
    params['fu']=fixed_update
    params['ri_smoothing']=ri_smoothing
    params['update_factor']=update_factor
    params['FW_tol']=FW_tol
    params['tol']=tol

    #Save the different documents
    to_save=[G_FW, OD, ri_FW, n_outer, n_inner, balance, opt_res, OD_list, balance_list, params]

    print('Saving to external .pkl file')
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    with open(os.path.join(save_dir,'output'+'_L_'+str(L_r)+'_ni_'+str(ni)+'_no_'+str(no)+'_ev_'+str(ev)+
    '_'+stopping_criterion+'_fu_'+str(fixed_update)+'.pkl'), 'wb') as f:
            pickle.dump(to_save, f)

    return


if __name__ == "__main__":
    main()
