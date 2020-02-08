import numpy as np 
import networkx as nx
import pandas as pd
import pickle
import argparse
import os
from result_analysis import print_balance, plot_ri
import matplotlib.pyplot as plt
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", help="path to data directory.",
                        type=str)
            
    parser.add_argument("-pc", help="path of files for comparison") 

    args = parser.parse_args()
    path = os.path.join('Data/',args.p,'outputs',args.pc)

    _ = plt.figure(figsize=(10,10))
    # ax1 = fig.add_subplot(1, 1, 1)
    for f in os.listdir(path):
        if not f.endswith('.pkl'):
            continue
        f_parsed=f.split(".")[0]
        f_parsed=f_parsed.split("_")
        L_r=f_parsed[2]
        ni=f_parsed[4]
        no=f_parsed[6]

        if args.pc.startswith("L"):
            value=L_r
        elif args.pc.startswith("ni"):
            value=ni
        elif args.pc.startswith("no"):
            value=no
        else:
            print("The folder does not correspond to an already defined graph to draw")
            return
        # Default values
        with open(os.path.join(path,f), 'rb') as filename:
            G_FW, OD, ri_FW, n_outer, n_inner, balance = pickle.load(filename)

        #Generate balance plot
        balance_norm=np.linalg.norm(balance,axis=1)
        plt.plot(np.arange(0,balance_norm.shape[0],1), balance_norm, label=value)
    plt.grid(True)
    plt.yscale('log')
    plt.xlabel('Iteration #')
    plt.legend()
    plt.savefig(os.path.join(path,'balance.png'))
    return

if __name__ == "__main__":
    main()