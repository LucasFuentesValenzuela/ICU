#all the elements necessary to build the graph 
#essentially to perform appropriate graph expansion
import os
import numpy as np
import pandas as pd
from helpers_icu import phi
import networkx as nx
from routines_icu import update_costs
from FW_OuterUpdate import init_flows

def construct_graph(path):

    edges, OD = expand_graph(path)
    G=init_graph(edges)
    G=initEdgeAttr(G,edges)
    G=initNodeAttr(G,edges)
    G=init_flows(G,OD) #check which version to consolidate
    G=update_costs(G) 

    return G, OD

def init_graph(edges):

    #extract all nodes from the file, after expansion
    nodes_list=np.unique(edges['head'].tolist()+edges['tail'].tolist())
    #extract all the edges from the file, after expansion
    edge_list=[]
    for i in range(edges.shape[0]):
        edge_list.append((edges.loc[i,'head'],edges.loc[i,'tail']))
    #initialize the graph
    G=nx.DiGraph()
    G.add_nodes_from(nodes_list)
    G.add_edges_from(edge_list)

    return G

def expand_graph(path):
    
    #read external Excel files containing the data
    edges=pd.read_excel(os.path.join(path,'edges.xlsx'))
    OD_xl=pd.read_excel(os.path.join(path,'OD.xlsx'))

    #convert into strings for the first two columns
    for i in range(edges.shape[0]):
        for col in ['head', 'tail']:
            edges.loc[i,col]=str(edges.loc[i,col])
    
    edges=edges.copy()
    
    o_list=np.unique(OD_xl['origin'])
    # d_list=np.unique(OD_xl['destination'])
    
    #######################################
    # 1 . Add rebalancing edge
    
    #what really matters is the value of phi
    L_rebalancing_edge=10
    t_rebalancing_edge=1
    k_rebalancing_edge=1 #does not matter as it should be adapted accordingly
    
    for o in o_list: #loop over the different origins in OD pairs
        data=np.array([[str(o),'R', L_rebalancing_edge,
                       k_rebalancing_edge,t_rebalancing_edge,0, 0]])
        l=pd.DataFrame(columns=edges.columns,data=data)
        edges=edges.append(l,ignore_index=True)
     
    
    #########################################################
    # 2. Add dummy destination nodes and inverse demand edges

    OD=dict()
    for i in range(OD_xl.shape[0]):
        #we add one dummy destination node per origin
        o=str(OD_xl.loc[i,'origin'])
        d=str(OD_xl.loc[i,'destination'])
        demand=OD_xl.loc[i,'demand']
        attr_=OD_xl.iloc[i,3:].values
        d_dummy=str(o)+'_p'
        
        if (o,d_dummy) in OD.keys():
            OD[o,d_dummy]+=demand
        else:
            OD[o,d_dummy]=demand
            
        data=np.array([[d,d_dummy]+attr_.tolist()])
        l=pd.DataFrame(columns=edges.columns,data=data)
        edges=edges.append(l,ignore_index=True)


    #######################################
    # 3. Add zero-cost edge

    #what really matters is the value of phi (0)
    L_ZC=0
    t_ZC=1
    k_ZC=1
    for o in o_list:
        data=np.array([[str(o), str(o)+'_p', L_ZC, k_ZC, t_ZC, 0, 0]])
        l=pd.DataFrame(columns=edges.columns,data=data)
        edges=edges.append(l,ignore_index=True)

    #make sure we have floats and no strings for the numbers
    edges.iloc[:,2:]=edges.iloc[:,2:].astype(float)
    return edges, OD#, dummy_nodes # I am not sure we actually need the dummy_node construct


def initEdgeAttr(G,edges):

    #initEdgeAttributes using a DataFrame

    for i in range(edges.shape[0]):
        e=(edges.loc[i,'head'],edges.loc[i,'tail'])
        G[e[0]][e[1]]['k']=edges.loc[i,'capacity']
        L=edges.loc[i,'length']
        t=edges.loc[i,'time']
        shift=edges.loc[i,'shift']
        G[e[0]][e[1]]['phi']=phi(L,t)
        is_negative=edges.loc[i,'isnegative']
        G[e[0]][e[1]]['sign']=(-1)**is_negative
        G[e[0]][e[1]]['f_m']=0
        G[e[0]][e[1]]['f_r']=0 
        G[e[0]][e[1]]['shift']=shift
    return G

def initNodeAttr(G,edges):

    max_shift=np.max(edges['shift'])

    for n in G.nodes():
        if n.endswith('_p'):
            G.nodes[n]['pot']=np.around(max_shift*1.1,0)
        G.nodes[n]['ri']=0

    return G

def vect_attribute(G,att):
    # x=[]
    x=dict()
    for e in G.edges():
        # x.append(G[e[0]][e[1]][att])
        x[e]=G[e[0]][e[1]][att]
    return x

def get_edge_list(G):
    l=[]
    for e in G.edges():
        l.append(e)

    return l

