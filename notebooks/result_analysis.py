import numpy as np
import matplotlib.pyplot as plt
from helpers_icu import BPR
import networkx as nx

def plot_edge_attrs(G_list, y_list, attrs, dots=True, lims=None):
    G_ = G_list[0]
    _, axes = plt.subplots(len(G_.edges()), len(
        attrs), figsize=(20, 5*len(G_.edges())))
    i = 0

    if dots:
        options = '--o'
    else:
        options = '--'

    for e in G_.edges():
        for j in range(len(attrs)):
            att = []
            if attrs[j] == "y_m":
                for y_k in y_list:
                    att.append(y_k[e, 'f_m'])
            elif attrs[j] == "y_r":
                for y_k in y_list:
                    att.append(y_k[e, 'f_r'])
            else:
                for G in G_list:
                    att.append(G[e[0]][e[1]][attrs[j]])
            axes[i, j].plot(att, options, label=attrs[j])
            axes[i, j].grid(True)
            axes[i, j].set_xlabel('Iteration #')
            axes[i, j].set_title(' Edge : ' + str(e))
            axes[i, j].legend()
            axes[i, j].set_xlim(lims)
        i += 1


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
        c.append(((G[o][d]['cost']+G[d][d+'_p']['cost']-G[o]
                   [d+'_p']['cost'])/G[o][d+'_p']['cost']))
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
        c.append(G[o][d]['cost']+G[d][d+'_p']['cost'])
        f_m.append(G[o][d]['f_m'])
        f_r.append(G[o][d]['f_r'])
    c = np.abs(np.array(c))

    #drop the first component
    c = c[1:]

    f_m = np.array(f_m)
    f_r = np.array(f_r)
    fig, ax1 = plt.subplots(figsize=(18, 5))
    ax1.plot(c, 'ro--', label="cost")
    ax1.plot(np.linspace(1, c.shape[0], 50), tgt_cost *
             np.ones(50), 'g--', label='Target value')
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
    ax1.set_ylim([88, 92])
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
    for (o, d) in OD.keys():
        x = np.linspace(0, OD[o, d], 100)
        cost_list = []
        cost_ID_list = []
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
        dict_cost[o, d] = cost_list
        dict_cost_ID[o, d] = cost_ID_list

    return dict_cost, dict_cost_ID


def plot_cost_all_path(G, OD, o, d):
    cost, cost_ID = get_cost_all_path(G, OD)

    if (o, d) not in cost.keys():
        print("Wrong set of o,d chosen")
        return

    nplots = len(cost[o, d])
    _, axes = plt.subplots(nplots, 1, figsize=(10, 4*nplots))
    max_c = 0
    for c in cost[o, d]:
        max_c = np.max([max_c, c[-1]])
    for i in range(nplots):
        crt_ = cost[o, d][i]
        crt_ID = cost_ID[o, d][i]
        diff = np.abs(np.array(crt_)-np.array(crt_ID))
        index_min = np.argmin(diff)
        x = np.linspace(0, OD[o, d], 100)
        axes[i].plot(x, crt_, label="cost")
        axes[i].plot(x, cost_ID[o, d][i], label="ID")
        axes[i].grid(True)
        axes[i].legend()
        axes[i].set_ylim([0, 1.1*max_c])
        axes[i].set_title("origin: " + o + ", destination: " + d + " | intersect at x = " +
                          str(np.around(x[index_min], 2)) + " | cost: " + str(crt_[index_min]))

    return


def plot_inv_demand(N, phi, k, shift):
    x = np.linspace(0, N, 100)
    INV = -BPR(phi, x, k)+shift
    plt.figure()
    plt.plot(x, INV)
    plt.grid(True)
    plt.ylim([0, INV[0]*1.2])
    return



###############################
#
# Flow conservation at nodes
#

def check_flow_cons(G, OD): 
    """ 
    Loop through the nodes and compute the total incoming and outgoing flow. 
    It should cancel at every origin/destination node if the network is properly rebalanced. 
    """
    net_flow=dict()
    for n in G.nodes():
       net_flow[n]=0 
    for e in G.edges():
        o=e[0]
        d=e[1]
        if d.endswith('_p') or d=='R':#we are on a dummy edge on the graph
            continue
        net_flow[o]-=G[o][d]['f_m']+G[o][d]['f_r']
        net_flow[d]+=G[o][d]['f_m']+G[o][d]['f_r']
    
    #the set of nodes for which the flows should indeed be balanced. 
    l=[]
    for (o,d) in OD.keys():
        l.append(o)
        if d.endswith('_p'):
            d=d.split('_')[0]
        l.append(d)
    l=np.unique(np.array(l))
    return net_flow, l

def check_flow_cons_at_OD_nodes(G,OD):
    """
    Compute the flow balance only at nodes that below to OD. 
    """
    net_flow,l=check_flow_cons(G,OD)
    balance=[]
    for n in l:
        balance.append(net_flow[n])
    return np.array(balance)