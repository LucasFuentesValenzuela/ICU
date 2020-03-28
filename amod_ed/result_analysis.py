import numpy as np
import matplotlib.pyplot as plt
from amod_ed.helpers_icu import BPR, BPR_int_val
import networkx as nx
import os

def plot_edge_attrs(G_list, y_list, attrs, dots=True, lims=None, ri_=None, only_R=False):
    G_ = G_list[0]
    _, axes = plt.subplots(len(G_.edges()), len(
        attrs), figsize=(20, 7*len(G_.edges())))
    i = 0

    tot_cost=np.zeros((len(G_.edges()), len(G_list)))
    tot_cost_2=np.zeros((len(G_.edges()), len(G_list)))
    if dots:
        options = '--o'
    else:
        options = '--'

    for e in G_.edges():
        if only_R and e[1]!='R':
            continue
        for j in range(len(attrs)):
            att = []
            att_2 =[]
            if attrs[j] == "y_m":
                for y_k in y_list:
                    att.append(y_k[e, 'f_m'])
            elif attrs[j] == "y_r":
                for y_k in y_list:
                    att.append(y_k[e, 'f_r'])

            #pretty much copy-pasted from the value cost function
            elif attrs[j] == "tot_cost": #we reconstruct total value value
                for G in G_list:
                    x_k_e_m = G[e[0]][e[1]]['f_m']
                    x_k_e_r = G[e[0]][e[1]]['f_r']
                    phi = G[e[0]][e[1]]['phi']
                    k = G[e[0]][e[1]]['k']
                    if k < 10**-5:  # you eliminate the edges that are considered non-usable
                        att.append(np.nan)
                        att_2.append(np.nan)
                        # print("skipping because nan")
                        # print(G[e[0]][e[1]]['cost'])
                        continue

                    # I am assuming there will be syntaxic problems there
                    F_E =  BPR_int_val(phi, x_k_e_m + x_k_e_r, k)

                    #this has to be included because it is directly included in the definition of the cost function
                    if G[e[0]][e[1]]['sign'] == (-1):  # we have a negative edge
                        F_E -= (x_k_e_m + x_k_e_r)*G[e[0]][e[1]]['shift']  # INVERSE_DEMAND_SHIFT

                    # not entirely sure this needs to be here
                    if 'pot' in G.nodes[e[1]]:
                        F_E+=G.nodes[e[1]]['pot']*(x_k_e_m + x_k_e_r)
                    att.append(F_E)

                    try:
                        att_2.append(G[e[0]][e[1]]['tot_cost'])
                    except KeyError:
                        att_2.append(np.nan)

                tot_cost[i,:] = np.array(att)
                tot_cost_2[i,:] = np.array(att_2)
            else:
                for G in G_list:
                    att.append(G[e[0]][e[1]][attrs[j]])
                    if attrs[j]=='k':
                        if e[1]=='R':
                            att_2.append(ri_[e[0]])
            axes[i, j].plot(att, options, label=attrs[j])
            if attrs[j] =="tot_cost":
                axes[i,j].plot(att_2,options, label='from loop')
            if attrs[j]=='k' and not ri_==None:
                axes[i,j].plot(att_2,options, label='ri')
            axes[i, j].grid(True)
            axes[i, j].set_xlabel('Iteration #')
            if j==0:
                axes[i, j].set_title(' Edge : ' + str(e))
            axes[i, j].legend()
            axes[i, j].set_xlim(lims)
        i += 1
    return tot_cost, tot_cost_2

def plot_node_attrs(G_list, attrs, lims=None):
    G_ = G_list[0]
    _, axes = plt.subplots(len(G_.nodes()), len(
        attrs), figsize=(18, 5*len(G_.nodes())))
    i = 0
    for n in G_.nodes():
        for j in range(len(attrs)):
            att = []
            for G in G_list:
                att.append(G.nodes[n][attrs[j]])

            if len(attrs) == 1:
                axes[i].plot(att, '--', label=attrs[j])
                axes[i].grid(True)
                axes[i].set_xlabel('Iteration #')
                axes[i].set_title(' node : ' + str(n))
                axes[i].set_xlim(lims)

            else:
                axes[i, j].plot(att, '--o', label=attrs[j])
                axes[i, j].grid(True)
                axes[i, j].set_xlabel('Iteration #')
                axes[i, j].set_title(' node : ' + str(n))
                axes[i, j].legend()
        i += 1

def plot_ri(ri_FW, G_FW, lims=None, compare=True, path=None):
    """

    compare is a boolean variable that aims to indicate whether or not 
    we want to plot the comparison alongside (i.e. the real value at the end of the iterations)
    """
    nplots=len(ri_FW[0].keys())
    ncols=3
    nrows=int(np.ceil(nplots/ncols))
    _, axes = plt.subplots(nrows,ncols, figsize=(13, 5*nrows))
    ri_iter=dict()
    i=0
    j=0

    #ri_iter is the ri that has been prescribed as an argument to the optimization 
    for n in ri_FW[0].keys():
        ri_iter[n]=[]

    for ri in ri_FW:
        for n in ri.keys():
            ri_iter[n].append(ri[n])
    
    #compute the net rebalancing flow at each node
    #those are the real net flows that occur on the graph
    if compare:
        real_ri=dict()
        G_ri=dict()
        for n in G_FW[0].nodes():
            real_ri[n]=np.zeros((len(ri_FW),))
            G_ri[n]=np.zeros((len(ri_FW),))
        for i in range(len(G_FW)):
            G=G_FW[i]
            for e in G.edges():
                if e[1]!='R' and not e[1].endswith('_p'):
                    real_ri[e[0]][i]-=G[e[0]][e[1]]['f_r']
                    real_ri[e[1]][i]+=G[e[0]][e[1]]['f_r']
            for n in G.nodes():
                G_ri[n][i]=G.nodes[n]['ri']
    i=0
    for n in ri_iter.keys():
        
        if j==ncols:
            j=0
            i+=1
        if i==nrows:
            i=0
        axes[i,j].plot(np.arange(1,len(ri_iter[n]),1),ri_iter[n][:-1], 'o-', label='Goal. ri_t')
        if compare:
            axes[i,j].plot(np.arange(1,len(ri_iter[n]),1),real_ri[n][1:], 'o--', label='net flows')
            axes[i,j].plot(np.arange(1,len(ri_iter[n]),1),real_ri[n][1:]-ri_iter[n][:-1], 'o--', label='diff')

            # axes[i,j].plot(G_ri[n],'o--', label='G_ri')
        axes[i,j].grid(True)
        axes[i,j].set_xlabel('Iteration #')
        axes[i,j].set_title(' node : ' + str(n))
        axes[i,j].set_xlim(lims)
        axes[i,j].set_xticks(np.arange(1,len(ri_iter[n]),1))
        axes[i,j].legend()
        # axes[i,j].set_ylim([np.min(ri_iter[n][-5:])*.95 , np.max(ri_iter[n][-5:])*1.05])
        j+=1
    
    if path!=None:
        plt.savefig(os.path.join(path,'ri.png'))
    
    plt.close()
    return

def print_balance(balance, path=None):

    balance_norm=np.linalg.norm(balance,axis=1)
    plt.figure(figsize=(10,10))
    plt.plot(balance_norm)
    plt.grid(True)
    plt.yscale('log')
    plt.xlabel('Iteration #')
    plt.xticks(np.arange(0,balance_norm.shape[0],1))

    if path!= None:
        plt.savefig(os.path.join(path,'balance.png'))

    plt.close()
    return


def plot_stop_and_cost(opt_res, scale='log'):
    """
    Plots both the value of the stopping criterion and the total cost
    """
    nplots=len(opt_res)
    ncols=2
    nrows=int(np.ceil(nplots/ncols))
    _, axes = plt.subplots(nrows, ncols, figsize=(20, 5*nrows))
    for n in range(nplots):
        i = int(np.floor(n / ncols))
        j=n % ncols
        axes[i,j].plot(opt_res[n]['stop'], label='stopping criterion', color='b')
        axes[i,j].set_title('Outer loop # '+ str(n))
        axes[i,j].set_yscale(scale)
    #     axes[i,j].set_xscale('log')
        axes[i,j].grid(True)
        axes[i,j].set_ylabel('stop crit')
        ax_2=axes[i,j].twinx()
        ax_2.plot(opt_res[n]['obj'], label='Total cost', color='g')
        ax_2.set_ylabel('total cost')
        ax_2.set_yscale(scale)
        ax_2.set_yticks([np.min(opt_res[n]['obj']), np.max(opt_res[n]['obj'])])
        # axes[i,j].set_ylim([10**2,10**8])
        # axes[i,j].set_xlim([0,100])
        axes[i,j].legend()
        ax_2.legend()


#####################################################################
# #
# for quals
def plot_stop_and_cost_output(opt_res, path, scale='log', nframes=4):
    """
    Plots both the value of the stopping criterion and the total cost
    Version that saves it to a path and only prints the first x frames
    """
    import matplotlib
    font = {'size'   : 15, 'weight': 'normal'}
    matplotlib.rc('font', **font)

    nplots=len(opt_res)
    ncols=2
    nrows=int(np.ceil(nframes/ncols))
    _, axes = plt.subplots(nrows, ncols, figsize=(25, 5*nrows))
    for n in range(nframes):
        i = int(np.floor(n / ncols))
        j=n % ncols
        axes[i,j].plot(opt_res[n]['stop'], label='Stopping criterion', color='b')
        axes[i,j].set_title('Outer loop # '+ str(n))
        axes[i,j].set_yscale(scale)
    #     axes[i,j].set_xscale('log')
        axes[i,j].grid(True)
        axes[i,j].set_ylabel('Stopping Criterion')
        ax_2=axes[i,j].twinx()
        ax_2.plot(opt_res[n]['obj'], label='Total cost', color='g')
        ax_2.set_ylabel('Total cost')
        ax_2.set_yscale(scale)
        # axes[i,j].set_ylim([10**2,10**8])
        # axes[i,j].set_xlim([0,100])
        axes[i,j].legend(loc=1)
        ax_2.legend(loc=2)
    plt.savefig(path,transparent=True, dpi=400)


def plot_balance_list_output(balance_list, path, nframes=4):
    """
    Plots both the value of the stopping criterion and the total cost
    """
    import matplotlib
    font = {'size'   : 15, 'weight': 'normal'}
    matplotlib.rc('font', **font)

    nplots=len(balance_list)
    ncols=2
    nrows=int(np.ceil(nframes/ncols))
    _, axes = plt.subplots(nrows, ncols, figsize=(25, 5*nrows))
    for n in range(nframes):
        i = int(np.floor(n / ncols))
        j=n % ncols
        r_p = []
        for k in range(len(balance_list[n])-1):
            r_p.append(np.abs(balance_list[n][k]-balance_list[n][k+1])/balance_list[n][k])

        axes[i,j].plot(balance_list[n], label='balance', color='b')
        axes[i,j].set_title('Outer loop # '+ str(n))
        # axes[i,j].set_yscale('log')
    #     axes[i,j].set_xscale('log')
        axes[i,j].grid(True)
        axes[i,j].set_ylabel('Balance norm')
        # axes[i,j].legend()
        axes[i,j].set_ylim([0, 10])
        # ax_2=axes[i,j].twinx()
        # ax_2.plot(r_p, 'g')
        # ax_2.set_yscale('log')
        # ax_2.set_ylabel('relative change')

    plt.savefig(path,transparent=True, dpi=400)

#####################################################################

def plot_balance_list(balance_list, b_scale='linear', progress = False):
    """
    Plots both the value of the stopping criterion and the total cost
    """
    nplots=len(balance_list)
    ncols=2
    nrows=int(np.ceil(nplots/ncols))
    _, axes = plt.subplots(nrows, ncols, figsize=(22, 5*nrows))
    for n in range(nplots):
        i = int(np.floor(n / ncols))
        j=n % ncols
        r_p = []
        for k in range(len(balance_list[n])-1):
            r_p.append(np.abs(balance_list[n][k]-balance_list[n][k+1])/balance_list[n][k])

        axes[i,j].plot(balance_list[n][1:], label='balance', color='b')#we put [1:] because want not to show drop in the beginning: TODO: understand fully and explain
        axes[i,j].set_title('Outer loop # '+ str(n))
        # axes[i,j].set_yscale('log')
    #     axes[i,j].set_xscale('log')
        axes[i,j].grid(True)
        axes[i,j].set_ylabel('balance norm')
        axes[i,j].legend()
        axes[i,j].set_yscale(b_scale)
        if progress: 
            ax_2=axes[i,j].twinx()
            ax_2.plot(r_p, 'g')
            ax_2.set_yscale('log')
            ax_2.set_ylabel('relative change')

def plot_ri_list(ri_FW, save = False, path = None):
    import matplotlib
    font = {'size'   : 15, 'weight': 'normal'}
    matplotlib.rc('font', **font)
    r=dict()
    for n in ri_FW[0].keys():
        r[n]=[]
        
    for n in r.keys():
        for ri in ri_FW:
            r[n].append(ri[n])

    plt.figure(figsize=(13,10))
    for n in r.keys():
        plt.plot(r[n],'--', label=n)
    plt.grid()
    # plt.legend()
    plt.xlabel('Iteration #')
    plt.ylabel('Rebalancing request per node')
    
    if save: 
        plt.savefig(path, transparent=True, dpi=400)


def compare_ri(ri_1, ri_2, node):
    r1=dict()
    r2=dict()
    for n in ri_1[0].keys():
        r1[n]=[]
        r2[n]=[]
        
    for n in r1.keys():
        for ri in ri_1:
            r1[n].append(ri[n])
        for ri in ri_2:
            r2[n].append(ri[n])

    plt.figure(figsize=(13,10))
    plt.plot(r1[node],'-', label='1')
    plt.plot(r2[node], '--', label='2')
    plt.grid()
    plt.legend()


def plot_OD(OD_list, o, d):
    vals = []
    for OD in OD_list:
        vals.append(OD[o, d])
    plt.figure(figsize=(13, 5))
    plt.plot(np.array(vals), 'o--', markersize=2)
    plt.grid(True)

####################################################################
########### DEBUG HELPERS
####################################################################

def disp_costs(G):

    print("Cost for dummy edge: ", G['1']['2_p']['cost'])
    print("Cost for normal edge (1,2),(2,2_p): ",
          G['1']['2']['cost'], " --- ", G['2']['2_p']['cost'])
    print("Cost for normal edge (1,2,2_p): ",
          G['1']['2']['cost']+G['2']['2_p']['cost'])
    print("Cost for 1-R: ", G['1']['R']['cost'])
    print("Cost for 2-R: ", G['2']['R']['cost'])


def print_final_flows(G_k):
    G_end = G_k[-1]
    for e in G_end.edges():
        print(e, " : ", G_end[e[0]][e[1]]['f_m']+G_end[e[0]][e[1]]['f_r'])


def print_final_cost(G_k):
    G_end = G_k[-1]
    for e in G_end.edges():
        print(e, " : ", G_end[e[0]][e[1]]['cost'])


def sanity_check_N(G_k):
    G_end = G_k[-1]
    print("1 to 2: ", G_end['1']['2_p']['f_m']+G_end['1']['2']['f_m'])
    print("2 to 1: ", G_end['2']['1_p']['f_m']+G_end['2']['1']['f_m'])


def sanity_check_cost(G_k):
    G_end = G_k[-1]
    print("1 to 2: ", G_end['1']['2_p']['cost'], " ===== ",
          G_end['1']['2']['cost']+G_end['2']['2_p']['cost'])
    print("2 to 1: ", G_end['2']['1_p']['cost'], " ===== ",
          G_end['2']['1']['cost']+G_end['1']['1_p']['cost'])


def plot_errors(G_k, scale='log', dots=True, lims=None):
    c_12 = []
    c_21 = []
    for G in G_k:
        c_12.append(np.abs((G['1']['2']['cost']+G['2']['2_p']
                            ['cost']-G['1']['2_p']['cost'])/G['1']['2_p']['cost']))
        c_21.append(np.abs((G['2']['1']['cost']+G['1']['1_p']
                            ['cost']-G['2']['1_p']['cost'])/G['2']['1_p']['cost']))

    if dots:
        options = '--o'
    else:
        options = '-'
    _, axes = plt.subplots(2, 1, figsize=(13, 10))
    axes[0].plot(c_12, options, label="(1,2)")
    axes[1].plot(c_21, options, label="(2,1)")
    for i in [0, 1]:
        axes[i].legend()
        axes[i].grid(True)
        axes[i].set_yscale(scale)
        axes[i].set_ylabel("relative error")
        axes[i].set_xlabel("Iteration #")
        axes[i].set_xlim(lims)


def analyze_cost_oscillations(G_k, o, d, lims=None, scale='log'):
    c = []
    f_m = []
    f_r = []
    for G in G_k:
        c.append(((G[o][d]['cost']+G[d][o+'_p']['cost']-G[o][o+'_p']['cost'])/G[o][o+'_p']['cost']))
        f_m.append(G[o][d]['f_m'])
        f_r.append(G[o][d]['f_r'])
    c = np.abs(np.array(c))
    f_m = np.array(f_m)
    f_r = np.array(f_r)
    fig, ax1 = plt.subplots(figsize=(18, 5))
    ax1.plot(c, 'ro', label="("+o+","+d+")")
    # ax1.plot(-c,'ro')
    ax2 = ax1.twinx()
    ax2.plot(f_m, "--o", label="f_m")
    ax2.plot(f_r, "--o", label="f_r")
    ax2.plot(f_m+f_r, "--o", label="sum")
    fig.legend()
    ax1.grid(True)
    ax1.set_yscale(scale)
    ax1.set_xlabel("Iteration #")
    ax1.set_ylabel("error")
    ax2.set_ylabel("Flow")
    if not lims == None:
        ax1.set_xlim(lims)


def analyze_cost_oscillations_2(G_k, o, d, lims=None, scale='linear'):
    tgt_cost = 90
    tgt_flow = 8.65
    c = []
    f_m = []
    f_r = []
    for G in G_k:
        c.append(G[o][d]['cost']+G[d][o+'_p']['cost'])
        f_m.append(G[o][d]['f_m'])
        f_r.append(G[o][d]['f_r'])
    c = np.abs(np.array(c))

    #drop the first component
    c = c[1:]

    f_m = np.array(f_m)
    f_r = np.array(f_r)
    fig, ax1 = plt.subplots(figsize=(8, 5))
    ax1.plot(c, 'ro--', label="cost")
    # ax1.plot(np.linspace(1, c.shape[0], 50), tgt_cost *
    #          np.ones(50), 'g--', label='Target value')
    # ax1.plot(-c,'ro')
    ax2 = ax1.twinx()
    # ax2.plot(f_m,"--o",label="f_m")
    # ax2.plot(f_r,"--o",label="f_r")
    ax2.plot(f_m+f_r, "--o", label="total flow (m + r)")
    fig.legend()
    ax1.grid(True)
    ax1.set_yscale(scale)
    ax1.set_xlabel("Iteration #")
    ax1.set_ylabel("Total cost")
    ax2.set_ylabel("Flow")
    tgt_flow=np.mean(f_m[-40:]+f_r[-40:])
    ylim_1=np.mean(c[-40:])
    ax1.set_ylim([ylim_1*0.99, ylim_1*1.01])
    ax2.set_ylim([tgt_flow*0.98, tgt_flow*1.02])

    if not lims == None:
        ax1.set_xlim(lims)


def analyze_cost_oscillations_3(G_k, y_k, o, d, lims=None, scale='linear'):
    #same analysis as above but focusing on the dummy edges

    tgt_cost = 90
    # tgt_flow=8.65
    c = []
    f_m = []
    f_r = []
    y_dummy = []
    for G in G_k:
        c.append(G[o][d]['cost']+G[d][d+'_p']['cost'])
        f_m.append(G[o][d+"_p"]['f_m'])
        f_r.append(G[o][d+"_p"]['f_r'])
    for y in y_k:
        y_dummy.append(y[(o, d+"_p"), 'f_m'])

    #drop the first component
    c = c[1:]

    c = np.abs(np.array(c))
    f_m = np.array(f_m)
    f_r = np.array(f_r)
    y_dummy = np.array(y_dummy)
    fig, ax1 = plt.subplots(figsize=(18, 5))
    # ax1.plot(-c,'ro')
    ax2 = ax1.twinx()
    ax1.plot(c, 'ro', label="Cost of the normal edge")
    ax1.plot(np.linspace(1, c.shape[0], 50), tgt_cost *
             np.ones(50), 'g--', label='Target value')
    ax2.plot(f_m, "--o", label="f_m dummy edge")
    ax2.plot(f_r, "--o", label="f_r dummy edge")
    for k in range(y_dummy.shape[0]):
        if y_dummy[k] != 0:
            ax1.plot(k*np.ones(50), np.linspace(1, 100, 50),
                     'k-', linewidth=1.3)
    # ax2.plot(y_dummy,label='y_m')
    fig.legend()
    ax1.grid(True)
    ax1.set_yscale(scale)
    ax1.set_xlabel("Iteration #")
    ax1.set_ylabel("Total cost")
    ax2.set_ylabel("flow")
    ax1.set_ylim([88, 92])
    avg_ = np.mean(f_m[-7:-1]+f_r[-7:-1])
    ax2.set_ylim([(avg_)*0.95, (avg_)*1.05])
    if not lims == None:
        ax1.set_xlim(lims)

    print("Black lines indicate that the assignment variable y_m is activated for the dummy edge")

###############################
#
# Cost of all the paths between o and d
# and how they vary with the flow
#   good to get a rough idea of what the solution should be
#


def get_cost_all_path(G, OD):

    dict_cost = dict()
    dict_cost_ID = dict()
    paths=dict()
    for (o, d) in OD.keys():
        x = np.linspace(0, OD[o, d], 100)
        cost_list = []
        cost_ID_list = []
        path_list=[]
        for path in nx.all_simple_paths(G, source=o, target=d):
            # print("o, d", o, d)
            # print("path: ", path)
            cost = np.zeros(x.shape)
            cost_ID = np.zeros(x.shape)  # inverse demand
            for i in range(len(path)-1):
                e_0 = path[i]
                e_1 = path[i+1]
                phi = G[e_0][e_1]['phi']
                k = G[e_0][e_1]['k']
                cost_crt = BPR(phi, x, k)
                # if e_1.endswith('_p') and e_0!=e_1.split('_')[0]: #dummy node and not zero cost edge
                if G[e_0][e_1]['sign'] == -1:
                    cost_crt -= G[e_0][e_1]['shift']
                    cost_ID -= cost_crt
                    # cost_ID+=G.nodes[e_1]['pot']
                else:
                    cost += cost_crt
            cost_list.append(cost)
            cost_ID_list.append(cost_ID)
            path_list.append(path)
        dict_cost[o, d] = cost_list
        dict_cost_ID[o, d] = cost_ID_list
        paths[o,d] = path_list

    return dict_cost, dict_cost_ID, paths


def plot_cost_all_path(G, OD, o, d):
    cost, cost_ID, paths = get_cost_all_path(G, OD)

    if (o, d) not in cost.keys():
        print("Wrong set of o,d chosen")
        return

    nplots = len(cost[o, d])
    
    ncols=2
    nrows=int(np.ceil(nplots/2))
    _, axes = plt.subplots(nrows,ncols, figsize=(13, 4*nrows))
    max_c = 0
    for c in cost[o, d]:
        max_c = np.max([max_c, c[-1]])
    for i in range(nrows):
        for j in range(2):
            n=ncols*i+j
            if n >= len(cost[o,d]):
                break
            crt_ = cost[o, d][n]
            crt_ID = cost_ID[o, d][n]
            diff = np.abs(np.array(crt_)-np.array(crt_ID))
            index_min = np.argmin(diff)
            x = np.linspace(0, OD[o, d], 100)
            axes[i,j].plot(x, crt_, label="cost")
            axes[i,j].plot(x, cost_ID[o, d][n], label="ID")
            axes[i,j].grid(True)
            axes[i,j].legend()
            axes[i,j].set_ylim([0, 1.1*max_c])
            axes[i,j].set_title(str(paths[o,d][n])+ " | intersect at x = " +
                            str(np.around(x[index_min], 2)) + " | cost: " + str(np.around(crt_[index_min],2)))

    return


def plot_inv_demand(N, phi, k, shift):
    x = np.linspace(0, N, 100)
    INV = -BPR(phi, x, k)+shift
    plt.figure()
    plt.plot(x, INV)
    plt.grid(True)
    plt.ylim([0, INV[0]*1.2])
    return


