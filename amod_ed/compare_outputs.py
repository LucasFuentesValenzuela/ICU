import numpy as np


def compare_final(G_icu_list,G_FW_list,att):
    G_icu=G_icu_list[-1]
    G_FW=G_FW_list[-1]
    
    for e in G_icu.edges():
        if att=='flow':
            att_icu=G_icu[e[0]][e[1]]['f_m']+G_icu[e[0]][e[1]]['f_r']
            att_FW=G_FW[e[0]][e[1]]['f_m']+G_FW[e[0]][e[1]]['f_r']
        elif att=='cost':
            att_icu=G_icu[e[0]][e[1]]['cost']
            att_FW=G_FW[e[0]][e[1]]['cost']
        try:
            diff_rel=abs((att_FW-att_icu)/att_FW)
        except:
            if att_FW==0:
                if att_FW-att_icu==0:
                    diff_rel=0.0
            else:
                diff_rel=np.nan
        print(e," : ", np.around(100*diff_rel,2),"%")