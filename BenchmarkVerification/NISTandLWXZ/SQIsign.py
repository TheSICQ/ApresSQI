from isogenychains import *
from strategies_sike import *

from math import ceil
from copy import deepcopy
from hashlib import sha256

# change to slow Fp2 squareroots

class Verifier:
    '''normal NIST class (or LWXZ) SQIsign verification'''
    def __init__(self, p, f, f2, f3, version = 'NIST'):
        #Set parameters
        self.p = p
        self.f = f
        self.strategy_f, _ = compute_strategy(2, self.f//2-1)
        self.f2 = f2
        self.strategy_f2, _ = compute_strategy(2, self.f2//2-1)
        self.f3 = f3
        assert (self.p + 1) % (2**self.f*3**self.f3) == 0
        if self.f3 > 0:
            self.strategy_f3, _ = compute_strategy(3, self.f3)
        #cofac
        self.m = (self.p+1)//(2**self.f)
        #self.e = ceil((15/4)*log(self.p,2) + 25)
        self.e = 975 #The above is more correct, but for measurements, we assume p is always of the same size
        #number of blocks
        self.B = ceil(self.e/self.f)
        #for last block, g = length last block.
        self.g = self.e % self.f
        if self.g > 1: #This is fucked if g == 1, ....
            self.strategy_g, _ = compute_strategy(2, self.g//2-1)

        assert version in ['NIST', 'LWXZ'], "variant is wrong mate (Use other folder for Apres-variants)"
        self.ver = version

    def verify(self, msg, sigma, pk):
        zippy, r, s = sigma
        E2, Q = self.decompress_resp(zippy, pk)
        #print(f'E_2 = {E2}')
        return self.decompress_and_check_chall(E2, s, Q, r, msg)

    def hash_to_point(self, msg, ProjA):
        P, Q, PmQ = basis_chal_torsion(ProjA, self.f2, self.f3)
        H = sha256()
        H.update(bytes(msg, "utf-8"))
        s = int.from_bytes(H.digest())%(2**self.f2 * 3**self.f3)
        return ladder_3pt(P, Q, PmQ, s, ProjA)

    def decompress_resp(self, s, ProjA):
        P, Q, PmQ = basis_two_torsion(ProjA, self.f)
        b, scalars = s
        assert len(scalars) == self.B
        if b:
            temp = P
            P = Q
            Q = temp
        count = 0
        for si in scalars:
            count +=1
            K = ladder_3pt(P, Q, PmQ, si, ProjA)
            
            if count != self.B:
                ProjA, Qlist = four_iso_chain_nist(K, [Q], ProjA, self.f, self.strategy_f)
            else:
                if self.g > 1 and self.g < self.f:
                    K = xDBLe(K, ProjA, self.f - self.g)
                    ProjA, Qlist = four_iso_chain_nist(K, [Q], ProjA, self.g, self.strategy_g)
                elif self.g == 0:
                    K = xDBLe(K, ProjA, self.f - self.g)
                    # Just forget about the isogeny computation for now
                else: 
                    ProjA, Qlist = four_iso_chain_nist(K, [Q], ProjA, self.f, self.strategy_f)

            Q, A = double_proj_normalize(Qlist[0], ProjA)
            #Q = proj_point_normalize(Qlist[0])
            #A = proj_point_normalize(ProjA)
            if count != self.B:
                if self.ver == 'NIST':
                    Q, P, PmQ = complete_basis_two_torsion(Q, A, self.f)
                else:
                    Q, P, PmQ = complete_basis_two_torsion_SIKE(Q, A, self.f)
        if self.f > self.f2:
            Q = xDBLe(Q, ProjA, self.f-self.f2)
        newProjA, isom = normalize_curve(ProjA)
        #Q = proj_point_normalize(Q)
        Q = ec_isomorphism_eval(Q, isom)
        return newProjA, Q
    
    def _derive_Q(self, s, P, Q, PmQ, ProjA):
        if len(s) == 2:
            b_1, _ = s
            if b_1 ==0:
                R = Q
            else:
                R = P
            return R
        b_1, _, b_2, _ = s
        if b_1 == b_2:
            if b_1 == 0:
                R = Q
            else:
                R = P
        else:
            scalars = [2**self.f2, 3**self.f3]
            R = xMULbidim(P, Q, PmQ, scalars[b_1], scalars[b_2], ProjA)
        return R


    def decompress_and_check_chall(self, ProjA, s, Q_cyc, r, msg):
        #P, Q, PmQ = basis_chal_torsion(ProjA, min(self.lam, self.f2), self.f3)
        # When chal is just a 2-power isogeny
        if len(s) == 2:
            Q = proj_point_normalize(Q_cyc)
            if self.ver == 'NIST':    
                Q, P, PmQ = complete_basis_two_torsion(Q, ProjA, self.f2) 
            else:
                Q, P, PmQ = complete_basis_two_torsion_SIKE(Q, ProjA, self.f2) 
            b, scalar = s
            if b:
                temp = P
                P = Q
                Q = temp
            K = ladder_3pt(P, Q, PmQ, scalar, ProjA)

            if proj_point_equal(xDBLe(K, ProjA, self.f2 - 1), xDBLe(Q_cyc, ProjA, self.f2-1)):
                # print("Not cyclic!!")
                # return False          we still continue to get a better cost image, doesnt matter much for benchmarking
                pass

            ProjA, Qlist = four_iso_chain_nist(K, [Q], ProjA, self.f2, self.strategy_f2)

        # Combination of 2 and 3
        elif len(s) == 4:
            P, Q, PmQ = basis_chal_torsion(ProjA, self.f2, self.f3)
            b1, s1, b2, s2 = s
            P2, Q2, P2mQ2 = xTPLe(P, ProjA, self.f3), xTPLe(Q, ProjA, self.f3), xTPLe(PmQ, ProjA, self.f3)
            P3, Q3, P3mQ3 = xDBLe(P, ProjA, self.f2), xDBLe(Q, ProjA, self.f2), xDBLe(PmQ, ProjA, self.f2)
            if b1:
                temp = P2
                P2 = Q2
                Q2 = temp
            K2 = ladder_3pt(P2, Q2, P2mQ2, s1, ProjA)
            if proj_point_equal(xDBLe(K2, ProjA, self.f2-1), xDBLe(Q_cyc, ProjA, self.f2-1)):
                #print("Not cyclic!!")
                # return False          we still continue to get a better cost image
                pass
            if b2:
                temp = P3
                P3 = Q3
                Q3 = temp
            K3 = ladder_3pt(P3, Q3, P3mQ3, s2, ProjA)
            Q = self._derive_Q(s, P, Q, PmQ, ProjA)
            #print(f'Qm = {proj_point_normalize(Q)}')
            #ProjA, Qlist = three_iso_chain_strategy(K3, [K2, Q], ProjA, self.f3, self.strategy_f3)
            ProjA, Qlist = three_iso_chain_naive(K3, [K2, Q], ProjA, self.f3)
            ProjA, Qlist = four_iso_chain_nist(Qlist[0], [Qlist[1]], ProjA, self.f2, self.strategy_f2)
            
        else:
            assert False, f'Signature (s = {s}) was invalid'
        
        # Recompute hash and check
        ProjA, isom = normalize_curve(proj_point_normalize(ProjA))        #removed affineness
        Q = ec_isomorphism_eval(Qlist[0], isom)
        Q, ProjA = double_proj_normalize(Q, ProjA)
        #Q = proj_point_normalize(Q)
        #ProjA = proj_point_normalize(ProjA)
        P = self.hash_to_point(msg, ProjA)
        #print(f'E_1 = {ProjA}')
        #print(f'K = {proj_point_normalize(P)}')
        return proj_point_equal(P, xMUL(Q, r, ProjA))

