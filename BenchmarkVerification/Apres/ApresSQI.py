from isogenychains import *
from strategies_sike import *
from ec_misc import j_invariant

from math import ceil, log
from copy import deepcopy
from hashlib import sha256

class ApresSQI_verifier:
    '''ApresSQI class SQIsign verification'''
    def __init__(self, p, f, f2, f3):
        #D_chall = 2**f2*3**f3
        self.p = p
        self.f = f
        self.strategy_f, _ = compute_strategy(2, self.f//2)
        self.f2 = f2
        self.strategy_f2, _ = compute_strategy(2, self.f2//2)
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
        if self.g > 0:
            self.strategy_g, _ = compute_strategy(2, self.g//2)

    def verify(self, msg, sigma, pk, push=True):
        if len(sigma) == 3: #unseeded
            zip, r, s = sigma
            seeds = None
            seeds_ls = None
        elif len(sigma) == 4: #uncompressed
            return self.verify_uncompressed(msg, sigma, pk)
        elif len(sigma) == 5: #seeded
            zip, r, s, seeds, seeds_ls = sigma
        else:
            assert False, "Signature is corrupted"
        E_2, Q = self.decompress_resp(zip, pk, seeds = seeds, push = push)
        return self.decompress_and_check_chall(E_2, s, r, msg, seeds = seeds_ls, Q_push = Q)
    
    def verify_uncompressed(self, msg, sigma, pk):
        gens, E_1coeff, hash_seeds, _ = sigma
        E_1 = [E_1coeff, [1,0]]
        gens = [[Kx, [1,0]] for Kx in gens]
        A = pk
        for count, K in enumerate(gens):
            if count != self.B - 1 or self.g == 0:
                ProjA, _ = four_iso_chain_opt(K, [], A, self.f, self.strategy_f)
            else:
                ProjA, _ = four_iso_chain_opt(K, [], A, self.g, self.strategy_g)
            A = proj_point_normalize(ProjA)
        E_2 = proj_point_normalize(ProjA)
        P = self.hash_to_point(msg, E_1, seeds = hash_seeds)
        P = proj_point_normalize(P)
        E_2m, _ = four_iso_chain_opt(P, [], E_1, self.f2, self.strategy_f2)
        return j_invariant(E_2m) == j_invariant(E_2)
    
    def decompress_resp(self, s, ProjA, seeds = None, push = True):
        b, scalars = s
        assert len(scalars) == self.B
        A = ProjA
        if push:
            if not seeds is None:
                Q = [small_square_with_minus_one_nonsquare(seeds[0][1]),[1,0]]
                seeds[0] = seeds[0][0]
            else:
                Q = sample_improved_Q_plus_plus(A)
            for count, si in enumerate(scalars):
                if seeds:
                    A, Q = self.step(A, b, si, count, seed = seeds[count], Q_push = Q)
                else:
                    A, Q = self.step(A, b, si, count, Q_push = Q)
                b = 0
            if A[1] != [1,0]:
                return proj_point_normalize(A), Q
            else:
                return A, Q
        else:
            for count, si in enumerate(scalars):
                if not seeds is None:
                    A = self.step(A, b, si, count, seed = seeds[count])
                else:
                    A = self.step(A, b, si, count)
                b = 0
            if A[1] != [1,0]:
                return proj_point_normalize(A), None
            else:
                return A, None
        
    def decompress_and_check_chall(self, ProjA, s, r, msg, seeds = None, Q_push = None):
        A = ProjA
        if Q_push:
            Q = Q_push
            if A[1] != [1,0] or Q[1] != [1,0]:
                A, Q = double_proj_normalize(A, Q)
            if not seeds is None: #with seeds, with push
                P = [small_nonsquare(seeds[0]), [1,0]]
            else: #no seeds, with push
                P = sample_secondary_P(A)
        else:
            if A[1] != [1,0]:
                A = proj_point_normalize(A)
            if not seeds is None: #with seeds, no push        
                i, j = seeds[0]
                P = [small_nonsquare(i), [1, 0]]
                Q = [small_square_with_minus_one_nonsquare(j), [1, 0]]
            else: #no seeds, no push
                Q = sample_improved_Q_plus_plus(A) 
                P = sample_secondary_P(A)
        #print(f"E_2 ver: {A}")
        PmQ = point_difference(P, Q, A)
        Qm = self._derive_Q(s, P, Q, PmQ, A)
        # When chal is just a 2-power isogeny
        if len(s) == 2:
            b, scalar = s
            K = ladder_3pt(P, Q, PmQ, scalar, A)
            m2 = (self.p+1)//(2**self.f2)
            K = xMULaffs(K, m2, A)
            Qm = xMULaffs(Qm, m2, A)

            ProjA, Qlist = four_iso_chain_opt(K, [Qm], A, self.f2, self.strategy_f2)

        # Combination of 2-, and 3-power
        elif len(s) == 4:  
            b1, s1, b2, s2 = s
            K2 = ladder_3pt(P, Q, PmQ, s1, A)
            if not b2:
                K3 = ladder_3pt(P, Q, PmQ, s2, A)
            else:
                K3 = ladder_3pt(Q, P, PmQ, s2, A)

            m2 = (self.p+1)//(2**self.f2)
            m3 = (self.p+1)//(3**self.f3)
            K2 = xMULaffs(K2, m2, A)
            K3 = xMULaffs(K3, m3, A)
            Qm = xMULaffs(Qm, m2, A)

            ProjA, Qlist = four_iso_chain_opt(K2, [K3, Qm], A, self.f2, self.strategy_f2)
            A = proj_point_normalize(ProjA)
            ProjA, Qlist = three_iso_chain_strategy(Qlist[0], [Qlist[1]], A, self.f3, self.strategy_f3)
        else:
            assert False, f'Signature (s = {s}) was invalid'

        # Recompute hash and check
        ProjA, isom = normalize_curve(ProjA)        #removed affineness
        #print(f"E_1 verified: {ProjA}")
        Q = ec_isomorphism_eval(Qlist[0], isom)
        Q, ProjA = double_proj_normalize(Q, ProjA)
        hash_seeds = None
        if seeds:
            hash_seeds = seeds[1]
        P = self.hash_to_point(msg, ProjA, seeds = hash_seeds)
        return proj_point_equal(P, xMULaff(Q, r, ProjA))
        

    def hash_to_point(self, msg, ProjA, seeds = None):
        A = ProjA
        if seeds:
            P = [small_nonsquare(seeds[0]), [1, 0]]
            Q = [small_square_with_minus_one_nonsquare(seeds[1]), [1, 0]]
            PmQ = point_difference(P, Q, A)
        else:
            P, Q, PmQ = basis_chal_torsion_SIKE(ProjA, self.f2, self.f3)
        H = sha256()
        H.update(bytes(msg, "utf-8"))
        s = int.from_bytes(H.digest(), 'big')%(2**self.f2 * 3**self.f3)
        K = ladder_3pt(P, Q, PmQ, s, A)
        m = (self.p + 1)//(2**self.f2*3**self.f3)
        K = xMULaffs(K, m, A)
        return K
    
    def step(self, ProjA, b, si, count, seed = None, Q_push = None):
        A = ProjA
        if Q_push:
            Q = Q_push
            if A[1] != [1,0] or Q[1] != [1,0]:
                A, Q = double_proj_normalize(A, Q)
            if not seed is None: #with seed, with push
                P = [small_nonsquare(seed), [1, 0]]
            else: #no seed, with push
                P = sample_secondary_P(A)
            P = proj_point_normalize(P)
        else:
            if A[1] != [1,0]:
                A = proj_point_normalize(A)
            if not seed is None: #with seed, no push
                P = [small_nonsquare(seed[0]), [1, 0]]
                Q = [small_square_with_minus_one_nonsquare(seed[1]), [1, 0]]
            else: #no seed, no push
                Q = sample_improved_Q_plus_plus(A) 
                P = sample_secondary_P(A) 

        #print(f"Q above {proj_point_normalize(xMULaffs(Q, self.m*2**(self.f-1), A))}")
        PmQ = point_difference(P, Q, A)
        if b:
            temp = deepcopy(P)
            P = deepcopy(Q)
            Q = temp
            b = 0
        #print(f"ver {count}: {proj_point_normalize(P), proj_point_normalize(Q)}")
        K = ladder_3pt(P, Q, PmQ, si, A)
        K = xMULaffs(K, self.m, A)
        if count != self.B - 1 or self.g == 0:
            if Q_push:
                ProjA, Qlist = four_iso_chain_opt(K, [Q], A, self.f, self.strategy_f)
                return ProjA, Qlist[0]
            else:
                ProjA, _ = four_iso_chain_opt(K, [], A, self.f, self.strategy_f)
                return ProjA
        else:
            K = xDBLe_aff(K, A, self.f - self.g)
            if Q_push:
                ProjA, Qlist = four_iso_chain_opt(K, [Q], A, self.g, self.strategy_g)
                return ProjA, Qlist[0]
            else:
                ProjA, Q = four_iso_chain_opt(K, [], A, self.g, self.strategy_g)
                return ProjA    

    def _derive_Q(self, s, P, Q, PmQ, ProjA):
        if len(s) == 2:
            b_1, _ = s
            if b_1 == 0:
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