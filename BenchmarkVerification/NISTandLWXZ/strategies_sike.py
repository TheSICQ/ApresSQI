####################################################################################################################
################        This file contains functions to compute optimal strategies for isogenies      ##############
################           Code taken from SIBC,  implements the algorithm from the SIKE spec         ##############
####################################################################################################################


def compute_strategy(ell, n):

    assert(ell == 2 or ell ==3)

    if ell == 2:
        p = 5633    #cost for 2 doublings         
        q = 5461   #cost for 4-isogeny evaluation
        if n == 0: # edge case, will be handled by odd isogeny
            return [], 0

    if ell == 3:
        p = 5322    #cost for tripling        
        q = 5282   #cost for 3-isogeny evaluation

    S = {1:[]}
    C = {1:0 }
    for i in range(2, n+1):
        b, cost = min(((b, C[i-b] + C[b] + b*p + (i-b)*q) for b in range(1,i)), key=lambda t: t[1])
        S[i] = [b] + S[i-b] + S[b]
        C[i] = cost
    return S[n], C[n]

if __name__=="__main__":
    Sn, Cn = compute_strategy(2, 124)
    print(Sn)
    print(Cn)