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
    parser.add_argument("-s", help="shift to apply to data for viz", type=int, default=0)
    args = parser.parse_args()
    plot_balance(args)
    plot_duality_gap(args)
    return 


def plot_balance(args):
    path = os.path.join('Data/',args.p,'outputs',args.pc)
    shift=args.s
    fig = plt.figure(figsize=(10,10))
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
            flag='L'
        elif args.pc.startswith("ni"):
            value=ni
            flag='ni'
        elif args.pc.startswith("no"):
            value=no
            flag='no'
        else:
            print("The folder does not correspond to an already defined graph to draw")
            return
        # Default values
        with open(os.path.join(path,f), 'rb') as filename:
            G_FW, OD, ri_FW, n_outer, n_inner, balance, opt_res = pickle.load(filename)

        #Generate balance plot
        balance_norm=np.linalg.norm(balance,axis=1)/np.sqrt(balance.shape[1])
        plt.plot(np.arange(0,balance_norm.shape[0]-shift,1), balance_norm[shift:], 'o--',label=value)
    plt.grid(True)
    plt.ylabel('Balance norm')
    plt.yscale('log')
    plt.xlabel('Iteration #')
    plt.xlim([0,15])
    ax=fig.axes[0]
    handles, labels = ax.get_legend_handles_labels()
    # sort both labels and handles by labels
    labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
    new_labels=[]
    for i in range(len(labels)):
        new_labels.append( flag + ': ' + str(labels[i]))
    ax.legend(handles, new_labels)
    # plt.legend()
    plt.savefig(os.path.join(path,'balance.png'))
    return

def plot_duality_gap(args):
    path = os.path.join('Data/',args.p,'outputs',args.pc)
    for f in os.listdir(path):
        if not f.endswith('.pkl'):
            continue
        with open(os.path.join(path,f), 'rb') as filename:
            _, _, _, _, _, _, opt_res = pickle.load(filename) 
        nplots=len(opt_res)
        ncols=3
        nrows=int(np.ceil(nplots/ncols))
        _, axes = plt.subplots(nrows, ncols, figsize=(20, 5*nrows))
        for n in range(nplots):
            i = int(np.floor(n / ncols))
            j=n % ncols
            axes[i,j].plot(opt_res[n]['dual_gap'])
            axes[i,j].set_title('Outer loop # '+ str(n))
            axes[i,j].set_yscale('log')
            axes[i,j].grid(True)
            axes[i,j].plot(opt_res[n]['obj'])
            axes[i,j].set_ylim([10**0,10**8])
        plt.savefig(os.path.join(path,f.split('.')[0]+'_dual_gap.png'))
    return

if __name__ == "__main__":
    main()