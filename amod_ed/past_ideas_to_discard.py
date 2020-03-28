

#Rebalancer smoothing that was in FW_outerupdate.


# some kind of rebalancer smoothing
# The idea is the following: if there is not enough progress, we believe
# it is because of the rebalancers
# therefore we start smoothing out only once there is not enough progress
# TODO: integrate in the estimate ri_k routine
r = dict()
for n in ri_k.keys():
    r[n] = []
for n in r.keys():
    for ri in ri_:
        r[n].append(ri[n])
        ind = np.minimum(1, len(r[n]))
        avg_n = np.mean(r[n][-ind:])
        # very strong!!
        if np.linalg.norm(balance_new) > np.linalg.norm(balance_k) and smoothing == False and i > 5:
            lim_i = i
            # smoothing=True
        if smoothing == True:
            # beta=2/(i+1-lim_i+2)
            beta = .5
            # print("ri smoothing on , beta: ", beta)
        else:
            beta = 1
        ri_k[n] = (beta)*ri_k[n]+(1-beta)*avg_n