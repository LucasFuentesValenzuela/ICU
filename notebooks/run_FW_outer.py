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
    parser.add_argument("path", help="path to data directory",
                        type=str)
            
    # parser.add_argument("path_out", help="path to output data directory",
    #                     type=str)

    args = parser.parse_args()
    path = os.path.join('Data/',args.path)

    if not os.path.exists(path):
        print("Path does not exist")
        return

    #TODO: include a possibility to tune the phi of rebalancing edge as a parameter
    #because we basically want to show the behavior with different values of phi
    G_0, OD = construct_graph(path)
    edge_list = get_edge_list(G_0)
    G_FW, ri_FW, n_outer, n_inner, balance = solve(
        G_0.copy(), OD.copy(), edge_list, tol=10**-2, FW_tol=10**-2, max_iter=5000)
    
    plot_ri(ri_FW, G_FW, lims=None, compare=True, path=path)
    print_balance(balance,path=path) 

    #Save the different documents
    to_save=[G_FW, OD, ri_FW, n_outer, n_inner, balance] 

    print('Saving to external .pkl file')
    with open(os.path.join(path,'output.pkl'), 'wb') as f:
            pickle.dump(to_save, f)

    return


if __name__ == "__main__":
    main()
