import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from graph import expand_graph
import pandas as pd
import os
from graph import construct_graph, get_edge_list
from FW_OuterUpdate import solve
from result_analysis import check_flow_cons_at_OD_nodes, check_flow_cons
import pickle
import argparse

def main():
    """
    Runs the algorithm on the datafiles given as arguments
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("path_in", help="path to data directory",
                        type=str)
            
    parser.add_argument("path_out", help="path to output data directory",
                        type=str)

    args = parser.parse_args()
    path = args.path_in

    G_0, OD = construct_graph(path)
    edge_list = get_edge_list(G_0)
    G_FW, ri_FW, n_outer, n_inner, balance = solve(
        G_0.copy(), OD.copy(), edge_list, tol=10**-4, FW_tol=10**-4, max_iter=3000)

    to_save=[G_FW, ri_FW, n_outer, n_inner, balance] 

    print('Saving to external .pkl file')
    with open(os.path.join(args.path_out,'output.pkl'), 'wb') as f:
            pickle.dump(to_save, f)

    return


if __name__ == "__main__":
    main()
