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

    args = parser.parse_args()
    path = os.path.join('Data/',args.p)
    save_dir=os.path.join(path,'outputs',args.sd)
    L_r=args.L
    ni=args.ni
    no=args.no
    ev=args.ev

    if ev==1:
        evolving_bounds=True
    else:
        evolving_bounds=False

    if not os.path.exists(path):
        print("Path does not exist")
        return


    print("---------------------------------------------------")
    print("Solving with the following parameters:")
    print("     path: ", path)
    print("     L_r: ", L_r)
    print("     ni: ", ni)
    print("     no: ", no)
    print("     ev: ", evolving_bounds)

    G_0, OD = construct_graph(path, L_rebalancing_edge = L_r)
    edge_list = get_edge_list(G_0)
    G_FW, ri_FW, n_outer, n_inner, balance, opt_res, OD_list = solve(
        G_0.copy(), OD.copy(), edge_list, tol=10**-2, FW_tol=10**-8, 
        max_iter_outer= no, max_iter= ni, evolving_bounds=evolving_bounds)
    
    #Save the different documents
    to_save=[G_FW, OD, ri_FW, n_outer, n_inner, balance, opt_res, OD_list]

    print('Saving to external .pkl file')
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    with open(os.path.join(save_dir,'output'+'_L_'+str(L_r)+'_ni_'+str(ni)+'_no_'+str(no)+'_ev_'+str(ev)+'.pkl'), 'wb') as f:
            pickle.dump(to_save, f)

    return


if __name__ == "__main__":
    main()
