import os
from graph import construct_graph, get_edge_list
from FW_OuterUpdate import solve
import pickle
import argparse

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

    G_0, OD = construct_graph(path)
    edge_list = get_edge_list(G_0)
    G_FW, ri_FW, n_outer, n_inner, balance = solve(
        G_0.copy(), OD.copy(), edge_list, tol=10**-4, FW_tol=10**-4, max_iter=3000)

    to_save=[G_FW, ri_FW, n_outer, n_inner, balance] 

    print('Saving to external .pkl file')
    with open(os.path.join(path,'output.pkl'), 'wb') as f:
            pickle.dump(to_save, f)

    return


if __name__ == "__main__":
    main()
