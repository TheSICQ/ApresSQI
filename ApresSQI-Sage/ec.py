from sage.all import *
from xonly import isMontgomery, xPoint, sqrtDeterministic

#############################
#                           #
#        Basis-stuff        #
#                           #
#############################
    
def CompleteBasis(R, D, x = 1):
    E = R.curve()
    i = E.base_field().gens()[0]
    cof = sqrt(E.order())/D 
    facD = factor(D)
    Drad = radical(D)
    Rsmalls = []
    Rsmall = R*(D/Drad)
    for (l, e) in facD:
        Rsmalli = Rsmall*(Drad/l)
        assert Rsmalli
        Rsmalls.append(Rsmalli)

    for _ in range(1000):
        x += i
        try:
            S = E.lift_x(x)*cof
        except:
            continue
        Ssmall = S*(D/Drad)
        basis = True
        for index, (l, e) in enumerate(facD):
            Ssmalli = Ssmall*(Drad/l)
            if (not Ssmalli) or (Rsmalls[index].weil_pairing(Ssmalli, l) == 1):
                basis = False
                break
        if basis:
            RmS = point_difference(R, S)
            if RmS != (R - S).xy()[0]:
                S = -S
                assert RmS == (R - S).xy()[0]
            S.set_order(D)
            return S
    assert False, "Something went wrong in Complete Basis..."

def TorsionBasis(E, D, xOnly = 0, seeds = None, small_ns = None, small_s = None):
    i = E.base_field().gens()[0]
    x = Integer(1)
    cof = sqrt(E.order())/D
    facD = factor(D)
    Drad = radical(D)    
    ## case 1: With seeds ##
    if seeds:
        p, q = seeds
        # Multiply by cofactor after, because of verification shenanigans
        P = E.lift_x(p)
        Q = E.lift_x(q)
        PmQ = point_difference(P, Q)
        if PmQ != (P - Q).xy()[0]:
            Q = -Q
            assert PmQ == (P - Q).xy()[0]
        P, Q = P*cof, Q*cof
        return P, Q
    ## case 2: Generate basis + seeds ##
    if small_ns and small_s:
        #P point
        for n, x in enumerate(small_ns):
            try:
                P = E.lift_x(x)
            except:
                continue
            Pexp = P*cof
            Psmall = Pexp*(D/Drad)
            fullorder = True
            for (l, e) in facD:
                if not Psmall*(Drad/l):
                    fullorder = False
                    break
            if fullorder:
                Pexp.set_order(D)
                break
        Psmalls = []
        for (l, e) in facD:
            Psmalli = Psmall*(Drad/l)
            assert Psmalli
            Psmalls.append(Psmalli)
        # Q point
        for m, x in enumerate(small_s):
            try:
                Q = E.lift_x(x)
            except:
                continue
            Qexp = Q*cof
            Qsmall = Qexp*(D/Drad)
            basis = True
            for index, (l, e) in enumerate(facD):
                Qsmalli = Qsmall*(Drad/l)
                if (not Qsmalli) or (Psmalls[index].weil_pairing(Qsmalli, l) == 1):
                    basis = False
                    break
            if basis:
                Qexp.set_order(D)
                break
        PmQ = point_difference(P, Q)
        if PmQ != (P - Q).xy()[0]:
            Qexp = -Qexp
        return Pexp, Qexp, (n, m)
    
    ## Case 3: No seeds involved
    while True:
        x += i
        try:
            P = E.lift_x(x)*cof
        except:
            continue
        Psmall = P*(D/Drad)
        fullorder = True
        for (l, e) in facD:
            if not Psmall*(Drad/l):
                fullorder = False
                break
        if fullorder:
            P.set_order(D)
            break
    Q = CompleteBasis(P, D)
    if xOnly == 1:
        return P, Q
    PmQ = P-Q
    xP, xQ, xPmQ = xPoint(P.xy()[0], E), xPoint(Q.xy()[0], E), xPoint(PmQ.xy()[0], E)
    if xOnly == 2:
        return P, Q, xP, xQ, xPmQ
    return xP, xQ, xPmQ


def point_difference(P, Q):
	#follows code from SQIsign NIST version 1.0
	#input: affine points P,Q
	#		affine curve constant ProjA = [A, [1,0]]
	#output: x-coordinate of P-Q = PmQ

	#check if all inputs are affine
    E = P.curve()
    assert isMontgomery(E)
    A = E.a2()
    Px, Qx = P.xy()[0], Q.xy()[0]
    PmQZ = Px - Qx
    t2 = Px*Qx
    t3 = t2 - 1
    t0 = PmQZ*t3
    PmQZ = PmQZ**2
    t0 = t0**2
    t1 = t2 + 1
    t3 = Px + Qx
    t1 = t1*t3
    t2 = t2*A
    t2 = 2*t2
    t1 = t1 + t2
    t2 = t1**2
    t0 = t2-t0
    assert t0.is_square()
    t0 = sqrtDeterministic(t0)
    PmQX = t0 + t1
    return PmQX/PmQZ

#############################
#                           #
#     Hashing to chal       #
#                           #
#############################

from hashlib import sha256

def hashToPoint(D, msg, E, seeds=None, small_ns=None, small_s=None):
    H = sha256()
    H.update(bytes(msg, "utf-8"))
    if small_ns and small_s:
        P, Q, seeds = TorsionBasis(E, D, small_ns=small_ns, small_s=small_s)
        s = int.from_bytes(H.digest())%D
        return P + s*Q, seeds
    else:
        if seeds:
            P, Q = TorsionBasis(E, D, seeds=seeds, small_ns=small_ns, small_s=small_s)
        else:
            P, Q = TorsionBasis(E, D, xOnly = 1)
        s = int.from_bytes(H.digest())%D
        return P + s*Q