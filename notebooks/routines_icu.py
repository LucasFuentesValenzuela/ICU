#routines to support in the numerical tests for the routing algorithm
from helpers_icu import BPR, BPR_int
import numpy as np

 


#sign not necessary (think about what is actually D-1 and about the fact that it has a negative sign in front)
def update_costs(G):
    UPPER_LIMIT=10**22

    for e in G.edges:
        x=G[e[0]][e[1]]['f_m']+G[e[0]][e[1]]['f_r']
        phi=G[e[0]][e[1]]['phi']
        k=G[e[0]][e[1]]['k']
        G[e[0]][e[1]]['cost']=BPR(phi,x,k)

        #if the capacity is too small (only happens for the edges to R)
        #then we can say that there is no flow there actually 
        # (there might have been flow allocated from previous iterations that will just
        # mess up our current estimations)
        if k<10**-5:#we eliminate the edges with a too high cost from the equation
            #TODO: figure out whether this is a problem for the iterative alg? 
            G[e[0]][e[1]]['cost']=UPPER_LIMIT
            G[e[0]][e[1]]['f_r']=0
            continue

        if G[e[0]][e[1]]['sign']==(-1): #we have a negative edge
            G[e[0]][e[1]]['cost']-=G[e[0]][e[1]]['shift']
            
        if 'pot' in G.nodes[e[1]]:
            G[e[0]][e[1]]['cost']+=G.nodes[e[1]]['pot']
    return G

def estimate_ri_k(G,ri_smoothing,a_k):
    #determine whether each node is in excess or deficit of rebalancers
    ri_k=dict()
    ri_k_prev=dict()
    
    if ri_smoothing:
        beta=a_k
    else:
        beta=1 #no smoothing

    for n in G.nodes():
        ri_k[n]=0
        ri_k_prev[n]=G.nodes[n]["ri"]
    for e in G.edges():
        if not e[1].endswith('_p') and e[1]!='R':
            ri_k[e[0]]+=G[e[0]][e[1]]['f_m']
            ri_k[e[1]]-=G[e[0]][e[1]]['f_m']
    for n in G.nodes():
        G.nodes[n]["ri"]= (1-beta) * ri_k_prev[n] + beta*ri_k[n]
    return ri_k,G

def update_OD(OD,ri_k, a_k, G, evolving_bounds=True):
    
    #update the OD pairs for rebalancers
    eps=10**-6
    for n in ri_k.keys():
        if n!='R' and not n.endswith('_p'):
            if ri_k[n]<-eps: #you are in excess
                OD[(n,'R')]=-ri_k[n]
            else:
                OD[(n,'R')]=0

    #update the bounds for "regular OD pairs"
    #todo: make sure that the costs are a good reflection of the flow at iteration k 
    #tocode: if x_k on inverse demand edge is too close to upper bound or lower bound you 
    #have to decrease or increase those bounds, respectively


    #TODO: improve the evolving bounds routine
    #currently, we keep lower bounds at zero
    if evolving_bounds:
        #parameters, to be given as arguments in the future
        #those parameters should really be given as arguments
        l1=10**-1
        l2=2 #such a high value currently disables it
        alpha=1.2 #I think a moving value on that would be better. And same for the l1/l2. should be based o2
        for (o,d) in OD.keys():
            
            if d.endswith('_p'): #only the nodes that are the dummy nodes
                crt_Ulim=OD[o,d]#currently, we treat only the upper limit 
                crt_p_flow = get_total_flow_to_dummy_node(G,d)
                rel_U_error=abs(crt_p_flow-crt_Ulim)/crt_p_flow
                # UPPER BOUND
                if rel_U_error <= l1:#x_k too close to U
                    new_Ulim=alpha*crt_Ulim 
                elif rel_U_error>=l2:#x_k too far from U
                    new_Ulim=alpha*crt_p_flow
                else:
                    new_Ulim=crt_Ulim
                OD[o,d]=new_Ulim 

    return OD

def update_capacities(G,ri_k): 
    eps=10**-6
    for n in G.nodes():
        if not n.endswith('_p') and n!='R':
            if ri_k[n]>eps: 
                G[n]['R']['k']=ri_k[n]
            else:
                G[n]['R']['k']=eps
    return G

#we currently assign a rebalancing flow to y_k, but this is useless currently
#as with the new cost function we do not take it into account
def AoN(G,OD):
    #perform the All Or Nothing assignment    
    y_k=init_y(G)
    eps=10**-6

    alpha=1.5 #TODO: what is this alpha for? ANSWER: for the update below (in comments)
    for (o,d) in OD.keys():
        U=OD[o,d]
        if U>eps:
            if d=='R':
                flag='f_r'
            else:
                flag='f_m'
            # print("AON, (o,d):", o,d)
            path=nx.shortest_path(G,source=o,target=d,weight='cost')

            #in this case, we are assigning to the dummy edge
            #again, not sure that the dummy edges are truly necessary in the end
            #TODO: figure out what the commented version here below actually meant
            #it seems to me like it was a fix to make sure we introduce the 
            #"complement" of the demand
            #we have to see whether or not it was helping in the ICU case... ? 
            if len(path)==2 and d!='R':
                # y_k[(o,d),'f_m']+=alpha*(U-G[dummy_nodes[d]][d]['f_m'])
                y_k[(o,d),'f_m']+=U
            else:
                for i in range(len(path)-1):
                    y_k[(path[i],path[i+1]),flag]+=U
    return y_k

def get_total_flow_to_dummy_node(G,d):
    """
    Computes the total flow going into a dummy node. 
    This is useful when updating the bound, as the bound is a measure of the max flow that can go into
    a dummy node. 
    """
    flow=0

    for e in G.edges():
        o_ = d.split('_')[0]
        if e[1]==d and e[0]!=o_:
            flow+=G[e[0]][d]['f_m']
    return flow

#TODO: consolidate in a single version the common functions to both FW outer and FW icu