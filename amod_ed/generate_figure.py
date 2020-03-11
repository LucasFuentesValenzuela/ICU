import numpy as np 
import networkx as nx
import pandas as pd
import pickle
import argparse
import os
from result_analysis import print_balance, plot_ri
import matplotlib.pyplot as plt



def main():
    import matplotlib
    font = {'size'   : 15, 'weight': 'normal'}
    matplotlib.rc('font', **font)

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", help="path to data directory.",
                        type=str)
            
    parser.add_argument("-pc", help="path of files for comparison") 
    parser.add_argument("-s", help="shift to apply to data for viz", type=int, default=0)
    args = parser.parse_args()
    path=args.p
    path_c=args.pc 
    shift=args.s
    plot_balance(path, path_c, shift)
    plot_duality_gap(path, path_c)
    return 


def plot_balance(path, path_c, shift, save=True):
    path = os.path.join('Data/', path,'outputs', path_c)
    fig = plt.figure(figsize=(10,8))
    # ax1 = fig.add_subplot(1, 1, 1)
    for f in os.listdir(path):
        print(f)
        if not f.endswith('.pkl'):
            continue
        f_parsed=f.split(".")[0]
        f_parsed=f_parsed.split("_")
        L_r=f_parsed[2]
        ni=f_parsed[4]
        no=f_parsed[6]

        if os.path.split(path_c)[-1].startswith("L"):
            value=L_r
            flag='L'
        elif os.path.split(path_c)[-1].startswith("ni"):
            value=ni
            flag='ni'
        elif os.path.split(path_c)[-1].startswith("no"):
            value=no
            flag='no'
        elif os.path.split(path_c)[-1].startswith("FW"):
            flag='FW'
        elif os.path.split(path_c)[-1].startswith("u"):
            value=no
            flag='no'
        else:
            print("The folder does not correspond to an already defined graph to draw")
            return
        # Default values
        with open(os.path.join(path,f), 'rb') as filename:
            G_FW, OD, ri_FW, n_outer, n_inner, balance, opt_res, OD_list, _, params = pickle.load(filename)
            # G_FW, OD, ri_FW, n_outer, n_inner, balance = pickle.load(filename)

        #Generate balance plot
        balance_norm=np.linalg.norm(balance,axis=1)/np.sqrt(balance.shape[1])
        if flag == 'FW':
            value=params['FW_tol']
        plt.plot(np.arange(0,balance_norm.shape[0]-shift,1), balance_norm[shift:], 'o--',label=value)
    plt.grid(True)
    plt.ylabel('Balance norm')
    plt.yscale('log')
    plt.xlabel('Iteration #')
    # plt.xlim([0,15])
    ax=fig.axes[0]
    handles, labels = ax.get_legend_handles_labels()
    # sort both labels and handles by labels
    labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
    new_labels=[]
    for i in range(len(labels)):
        new_labels.append( flag + ': ' + str(labels[i]))
    ax.legend(handles, new_labels)
    # plt.show()
    # plt.legend(handles, new_labels)
    if save:
        plt.savefig(os.path.join(path,'balance.png'),transparent=True, dpi=400)
    return

def plot_duality_gap(path, path_c):
    path = os.path.join('Data/',path,'outputs',path_c)
    for f in os.listdir(path):
        if not f.endswith('.pkl'):
            continue
        with open(os.path.join(path,f), 'rb') as filename:
            _, _, _, _, _, _, opt_res, _ = pickle.load(filename) 
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