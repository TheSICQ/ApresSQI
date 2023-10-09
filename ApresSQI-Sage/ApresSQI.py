from klpt import *
from id2iso import *
from param_loader import load, AlternativeOrder
from quaternion import pushforward, MakeCyclic
from ec import hashToPoint, CompleteBasis, TorsionBasis
from ast import literal_eval
from xonly import Normalized, xPoint

class SQIsign:
    """
    Class for the SQIsign digital signature scheme,
    using multiple field extensions when signing.
    """
    def __init__(self, param = 'toy'):
        params = load(param)
        self.p = params['p']
        self.F = params['F']
        self.f = params['f']
        self.D_com = params['D_com']
        self.D_chall = params['D_chall']
        self.T = params['T']
        self.B_2 = params['B_2']
        self.B_chall = params['B_chall']
        self.facToExt = params['facToExt']
        self.facToBasis = params['facToBasis']
        self.facToAction = params['facToAction']
        self.small_ns = params['small_ns']
        self.small_s = params['small_s']

        self.B = QuaternionAlgebra(-1, -self.p)
        i, j, k = self.B.gens()
        self.O0 = self.B.quaternion_order([1, i, (i+j)/2, (1+k)/2])
        self.O0.QF = BinaryQF([1,0,1])

        self.O0_alt = AlternativeOrder(self.p)

        self.E0 = EllipticCurve(self.F, [1,0])

        self.sk = None
        self.pk = None

    def KeyGen(self):
        print("Generating secret key ....")
        p_bitsize = ceil(log(self.p,2))
        I_secret_bitsize = ceil(p_bitsize/4) 
        # The secret degree
        D_secret = next_prime(randint(2**I_secret_bitsize, 2**(I_secret_bitsize+1)))
        while D_secret % 4 != 3:
            D_secret = next_prime(D_secret)
        # The secret ideal
        gamma = FullRepresentInteger(self.O0, D_secret * 2**p_bitsize)
        i,_,_ = self.B.gens()
        # The alternative path
        print("Starting KLPT")
        for _ in range(100):
            a = randint(0, D_secret)
            I_secret = self.O0*(gamma*(a + i)) + self.O0*D_secret
            alpha = KeyGenKLPT(self.O0, I_secret, self.f)
            if alpha:
                break
        print("KLPT done!")
        J_secret = I_secret*(alpha.conjugate()/D_secret)
        phi_secret, _, Q_2f = IdealToIsogenyEichler(self.O0, self.O0_alt, J_secret, self.O0*1, self.facToBasis, self.facToAction, self.B_2[0], self.f, self.T)
        E_A, phi_secret, Q_2f = Normalized(phi_secret, R = Q_2f)
        pushedFacToBasis = PushBasis(phi_secret, self.facToBasis)
        self.pk = SQIsign_pubkey(E_A)
        self.sk = SQIsign_privkey(alpha.conjugate(), pushedFacToBasis, Q_2f, self.pk)
        return self.pk

    def Sign(self, msg, seeded = True, compressed = True):
        B = self.O0.quaternion_algebra()
        i,j,k = B.gens()
        theta = j + (1+k)/2
        alpha_sec, pushedFacToBasis, Q_2f = self.sk.unpack()
        N_alpha_sec = alpha_sec.reduced_norm()
        pow2 = N_alpha_sec.valuation(2)
        D_secret = ZZ(N_alpha_sec/(2**pow2))
        J_secret = self.O0*alpha_sec + self.O0*(2**pow2)
        I_secret = self.O0*alpha_sec.conjugate() + self.O0*D_secret

        #<<<<---- Commitment --->>>>>
        #TODO: The commitment ideal is not correctly generated, figure how it should be (specs are completely wrong...)
        gamma = FullRepresentInteger(self.O0, self.D_com * 2**ceil(log(self.p, 2)))
        #a, b = [randint(1, self.D_com) for _ in range(2)]
        #I_com = self.O0*(a + b*theta) + self.O0*self.D_com
        I_com = self.O0*gamma + self.O0*self.D_com
        phi_com = IdealToIsogeny(self.O0, I_com, self.E0, self.facToBasis, self.facToAction)
        E_1, phi_com = Normalized(phi_com)
        P_1, Q_1 = [phi_com(P) for P in self.B_chall]

        #<<<<---- Challenge --->>>>>
        if seeded or not compressed:
            K_chall, seeds_hash = hashToPoint(self.D_chall, msg, E_1, small_ns = self.small_ns, small_s = self.small_s)
        else:
            K_chall = hashToPoint(self.D_chall, msg, E_1)
        phi_chall = E_1.isogeny(K_chall, algorithm='factored')
        E_2, phi_chall = Normalized(phi_chall)
        a, b = BiDLP(K_chall, P_1, Q_1, self.D_chall)
        I_chall = pushforward(I_com, KernelDecomposedToIdeal(self.O0, self.D_chall, a, b, self.facToAction))

        #<<<<---- Response --->>>>>
        K = I_secret.conjugate()*I_com*I_chall
        for _ in range(100):
            J = SigningKLPT(self.O0, K, I_secret, self.f)
            if not J:
                continue
            # test that composition with challenge is cyclic
            alpha_m = IdealEquivalence(J, K)
            J_temp = alpha_m**(-1)*J*alpha_m
            assert J_temp.right_order() == I_chall.right_order()
            if MakeCyclic(J_temp * I_chall.conjugate()) != J_temp * I_chall.conjugate():
                continue
            # Output ideal J
            alpha = IdealEquivalence(I_secret, J_secret)
            J = (alpha**(-1))*J*alpha
            # Only temporary, should not be necessary to enforce (actually, enforcing this leaks one secret bit lol):
            if MakeCyclic(J_secret*J) != J_secret*J:
                continue
            break
        assert J, "SigningKLPT failed"
        phi_sig, zip, _ = IdealToIsogenyEichler(self.O0, self.O0_alt, J, J_secret, pushedFacToBasis, self.facToAction, Q_2f, self.f, self.T, I_secret=I_secret)
        E_2 = phi_sig.codomain()

        #<<<<---- Recompute the hash --->>>>>
        seeds_last_step = []
        P = CompleteBasis(K_chall, self.D_chall)
        phi_chall = phi_chall.codomain().isomorphism_to(E_2) * phi_chall
        Pm = phi_chall(P)
        if seeded:
            P_2, Q_2, seeds_chall = TorsionBasis(E_2, self.D_chall, small_ns=self.small_ns, small_s = self.small_s)
            seeds_last_step.append(seeds_chall)
            seeds_last_step.append(seeds_hash)
        else:
            P_2, Q_2 = TorsionBasis(E_2, self.D_chall, xOnly=1)
        f2 = self.D_chall.valuation(2)
        f3 = self.D_chall.valuation(3)
        a_0, a_1 = BiDLP(Pm*3**f3, P_2*3**f3, Q_2*3**f3, 2**f2)
        b_1 = 0
        if a_0 % 2 == 0:
            temp = a_0
            a_0 = a_1
            a_1 = temp
            b_1 = 1
        s_1 = (pow(a_0, -1, 2**f2) * a_1) % 2**f2
        if f3 != 0:
            a_0, a_1 = BiDLP(Pm*2**f2, P_2*2**f2, Q_2*2**f2, 3**f3)
            b_2 = 0
            if a_0 % 3 == 0:
                temp = a_0
                a_0 = a_1
                a_1 = temp
                b_2 = 1
            s_2 = (pow(a_0, -1, 3**f3) * a_1) % 3**f3
            s = (b_1, s_1, b_2, s_2)
        else:
            s = (b_1, s_1)

        E_1m, phi_chall_hat = Normalized(E_2.isogeny(Pm, algorithm='factored'))
        Qm = self._derive_Q(s, P_2, Q_2)
        Q = phi_chall_hat(Qm)
        r = Q.discrete_log(K_chall)
        assert r*Q == K_chall
        if compressed:
            if seeded:
                newzip, seeds = self._getSeeds(self.pk.E_A, zip)
                sigma = (newzip, r, s, seeds, seeds_last_step)
            else:
                sigma = (zip, r, s)
        else:
            gens = [list(K.xy()[0]) for K in self._getGens(self.pk.E_A, zip)]
            sigma = (gens, list(E_1.a2()), seeds_hash, [])
        return sigma
    
    def verify(self, msg, sigma, pk):
        if len(sigma) == 3: #unseeded
            zip, r, s = sigma
            seeds = None
            seeds_ls = None
        elif len(sigma) == 4: #uncompressed
            return self.verify_uncompressed(msg, sigma, pk)
        elif len(sigma) == 5: #seeded
            zip, r, s, seeds, seeds_ls = sigma
        else:
            assert False, "Somethings wrong with the signature"
        E_2 = self.decompress_response(zip, pk.E_A, seeds = seeds)
        return self.decompress_and_check_chall(E_2, s, r, msg, seeds = seeds_ls)
    
    def decompress_response(self, s, E, seeds = None):
        if seeds:
            P, Q = TorsionBasis(E, 2**self.f, seeds = self._unpack_seeds(seeds[0]))
        else:
            P, Q = TorsionBasis(E, 2**self.f, xOnly = 1)
        b, scalars = s
        if b:
            temp = P
            P = Q
            Q = temp
        for i, si in enumerate(scalars):
            G = P + si*Q
            phi_i = customTwoIsogeny(G, self.f)
            if i < len(scalars) - 1:
                if seeds:
                    P, Q = TorsionBasis(phi_i.codomain(), 2**self.f, seeds = self._unpack_seeds(seeds[i+1]))
                else:
                    Q = phi_i(Q)
                    P = CompleteBasis(Q, 2**self.f)
        return phi_i.codomain()
    
    def decompress_and_check_chall(self, E, s, r, msg, seeds = None):
        if seeds:
            P_2, Q_2 = TorsionBasis(E, self.D_chall, seeds = self._unpack_seeds(seeds[0]))
        else:
            P_2, Q_2 = TorsionBasis(E, self.D_chall, xOnly = 1)        
        if len(s) == 2:
            b_1, s_1 = s
            f2 = Integer(self.D_chall).valuation(2)
            if b_1:
                Q_22 = P_2
                P_22 = Q_2
            else:
                P_22 = P_2
                Q_22 = Q_2
            K_22 = P_22 + s_1*Q_22
        else:
            b_1, s_1, b_2, s_2 = s
            f2 = Integer(self.D_chall).valuation(2)
            f3 = Integer(self.D_chall).valuation(3)
            if b_1:
                Q_22 = P_2*3**f3
                P_22 = Q_2*3**f3
            else:
                P_22 = P_2*3**f3
                Q_22 = Q_2*3**f3
            K_22 = P_22 + s_1*Q_22
            if b_2:
                Q_23 = P_2*2**f2
                P_23 = Q_2*2**f2
            else:
                P_23 = P_2*2**f2
                Q_23 = Q_2*2**f2
            K_23 = P_23 + s_2*Q_23
        Q = self._derive_Q(s, P_2, Q_2)
        phi_chall_hat = E.isogeny(K_22, algorithm='factored')
        if len(s) > 2:
            phi_chall_hat = phi_chall_hat.codomain().isogeny(phi_chall_hat(K_23), algorithm = 'factored')*phi_chall_hat
        E_1, phi_chall_hat = Normalized(phi_chall_hat)
        Q = phi_chall_hat(Q)
        if seeds:
            hash_seeds = self._unpack_seeds(seeds[1])
        else:
            hash_seeds = None
        K_chall = hashToPoint(self.D_chall, msg, E_1, seeds = hash_seeds)
        return K_chall.xy()[0] == (r*Q).xy()[0]
    
    def verify_uncompressed(self, msg, sigma, pk):
        gens, E_1coeff, hash_seeds, _ = sigma
        E_1 = EllipticCurve(self.F, [0, E_1coeff, 0, 1, 0])
        gens = [self.F(x) for x in gens]
        E_i = pk.E_A
        for xK in gens:
            K = E_i.lift_x(xK)
            phi_i = customTwoIsogeny(K, self.f)
            E_i = phi_i.codomain()
        E_2 = E_i
        K_chall = hashToPoint(self.D_chall, msg, E_1, seeds = self._unpack_seeds(hash_seeds))
        phi_chall = E_1.isogeny(K_chall, algorithm='factored')
        E_2m = phi_chall.codomain()
        return E_2m.j_invariant() == E_2.j_invariant()

    def _getGens(self, E_A, zip):
        r"""
        Recomputes the isogeny generators for the uncompressed version.
        Can of course also be done directly
        """
        gens = []
        E_i = E_A
        P, Q = TorsionBasis(E_i, 2**self.f, xOnly=1)
        b, scalars = zip
        if b:
            temp = P
            P = Q
            Q = temp
        G = P + scalars[0]*Q
        gens.append(G)
        phi_i = customTwoIsogeny(G, self.f)
        Q = phi_i(Q)
        P = CompleteBasis(Q, 2**self.f)
        for i, si in enumerate(scalars[1:]):
            G = P + si*Q
            E_i = phi_i.codomain()
            gens.append(G)
            phi_i = customTwoIsogeny(G, self.f)
            if i < len(scalars) - 2:
                Q = phi_i(Q)
                P = CompleteBasis(Q, 2**self.f)
        return gens
    
    def _derive_Q(self, s, P, Q):
        f2 = self.D_chall.valuation(2)
        f3 = self.D_chall.valuation(3)
        if len(s) == 2: #No 3 part
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
            scalars = [2**f2, 3**f3]
            R = scalars[b_1]*P + scalars[b_2]*Q
        return R
    
    def _getSeeds(self, E_A, zip):
        r"""
        Recompute the response, to get seeds for faster verification.
        Can be done directly, but for allowing non-seeded versions we do it after
        """
        newzip = []
        seeds = []
        E_i = E_A
        P, Q = TorsionBasis(E_i, 2**self.f, xOnly=1)
        b, scalars = zip
        if b:
            temp = P
            P = Q
            Q = temp
        G = P + scalars[0]*Q
        # First step a bit clunky
        newP, newQ, step_seeds = TorsionBasis(E_i, 2**self.f, small_ns = self.small_ns, small_s = self.small_s)
        seeds.append(step_seeds)
        a, b = BiDLP(G, newP, newQ, 2**self.f)
        new_b = 0
        if a % 2 == 0:
            new_b = 1
            temp = newP
            newP = newQ
            newQ = temp
            new_si = pow(b,-1,2**self.f)*a % 2**self.f
        else:
            new_si = pow(a,-1,2**self.f)*b % 2**self.f
        newzip.append(new_si)

        phi_i = customTwoIsogeny(G, self.f)
        Q = phi_i(Q)
        P = CompleteBasis(Q, 2**self.f)
        for i, si in enumerate(scalars[1:]):
            G = P + si*Q
            E_i = phi_i.codomain()
            newP, newQ, step_seeds = TorsionBasis(E_i, 2**self.f, small_ns = self.small_ns, small_s = self.small_s)
            seeds.append(step_seeds)
            a, b = BiDLP(G, newP, newQ, 2**self.f)
            new_si = pow(a,-1,2**self.f)*b % 2**self.f
            newzip.append(new_si)
            phi_i = customTwoIsogeny(G, self.f)
            if i < len(scalars) - 2:
                Q = phi_i(Q)
                P = CompleteBasis(Q, 2**self.f)
        return (new_b, newzip), seeds

    def _unpack_seeds(self, seeds):
        n, m = seeds
        small_ns = self.small_ns[n]
        small_s = self.small_s[m]
        return (small_ns, small_s)

    def load_privkey(self, param = None, fname = None):
        if param:
            filename = 'Privkeys/SQI_' + param + '_priv.key'
        else:
            assert fname, "Must give either file name or parameter name"
            filename = 'Privkeys/' + fname
        try:
            with open(filename, "r") as file:
                alpha_list = literal_eval(file.readline())
                pushedFacToBasis_str = literal_eval(file.readline())
                Q = literal_eval(file.readline())
                A = self.F(file.readline())
        except:
            assert False, 'no such file'
        d = alpha_list[-1]
        alpha = self.B(alpha_list[:-1])/d
        E_A = EllipticCurve(self.F, [0,A,0,1,0])
        Q_2f = E_A(Q)
        pushedFacToBasis = {}
        for fac in pushedFacToBasis_str.keys():
            Fbig = self.facToBasis[fac][0].curve.base_field()
            Ebig = E_A.base_extend(Fbig)
            P, Q, PmQ = [xPoint(Fbig(R), Ebig) for R in pushedFacToBasis_str[fac]]
            pushedFacToBasis[fac] = [P, Q, PmQ]
        
        self.pk = SQIsign_pubkey(E_A)
        self.sk = SQIsign_privkey(alpha, pushedFacToBasis, Q_2f, self.pk)

        return self.pk
    

class SQIsign_pubkey:
    """
    Class describing the SQIsign public key,
    which is just a montgomery coefficient.
    """
    def __init__(self, E_A):
        self.E_A = E_A


class SQIsign_privkey:
    """
    Class describing the SQIsign public key,
    which constist of an alpha describing the
    secret ideals, a pushed dictionary of
    torsion bases, a point Q useful in signing,
    and a SQIsign public key
    """
    def __init__(self, alpha, pushedFacToBasis, Q, pubkey):
        self.alpha = alpha
        self.pushedFacToBasis = pushedFacToBasis
        self.Q = Q
        self.pubkey = pubkey

    def unpack(self):
        return self.alpha, self.pushedFacToBasis, self.Q
    
    def export(self, fname):
        filename = 'Privkeys/' + fname
        with open(filename, 'w') as f:
            d = self.alpha.denominator()
            alpha_list = list((self.alpha*d).coefficient_tuple())
            alpha_list.append(d)
            f.write(f'{alpha_list}\n')
            printable_facToBasis = {}
            for key in self.pushedFacToBasis.keys():
                printable_facToBasis[key] = [[str(c) for c in P.X] for P in self.pushedFacToBasis[key]]
            f.write(f'{printable_facToBasis}\n') #f
            f.write(f'{[str(c) for c in self.Q.xy()]}\n')
            f.write(f'{self.pubkey.E_A.a2()}\n')