from sage.all import *
from quaternion import IdealEquivalence, generate_close_vectors, MakePrimitive, RandomEquivalentPrimeIdeal, IdealGenerator, QuaternionOrderBasis, pullback, ConnectingIdeal, IsomorphismGamma, EvalIsomorphism

def Cornacchia(QF, m):
    r"""
    Given a binaryQF QF, and a number m
    returns a solution of QF = m
    """
    m_prime = prod([l**e for l, e in factor(m, limit=100) if l < 100])
    if not is_pseudoprime(m/m_prime):
        return None, None, False
    sol = QF.solve_integer(m)
    if not sol:
        return None, None, False
    return sol[0], sol[1], True

def FullRepresentInteger(O0, M):
    r"""
    Given a spec. p-extremal order O0, and int M
    returns alpha \in O0 of norm M
    """
    B = O0.quaternion_algebra()
    QF = O0.QF
    i,j,k = B.gens()
    p = B.ramified_primes()[0]
    sol = None
    for _ in range(1000):
        m1 = max(floor(sqrt(round((4*M)/p, 5))), 100)
        z = randint(-m1, m1)
        m2 = floor(sqrt(round((4*M-z**2)/p, 5)))
        w = randint(-m2, m2)
        Mm = 4*M - p*QF(z,w)
        x, y, found = Cornacchia(QF, Mm)
        if not found:
            continue
        if x % 2 == w % 2 and y % 2 == z % 2 and gcd([x,y,z,w]) == 1:
            gamma = x + y*i + j*z + k*w
            gamma = MakePrimitive(gamma,O0)
            if gamma.reduced_norm() == M:
                return gamma
    return None


#####################################
#                                   #
#          Linear algebra           #
#                                   #
#####################################

def EichlerModConstraint(I, gamma, delta, divisible=False):
    r"""
    KLPT subroutine
    """
    B = I.quaternion_algebra()
    i,j,k = B.gens()
    beta1, beta2, _, _ = I.gens() # Should already be in HNF form
    if divisible:
        M = Matrix(Integers(I.norm()), [(gamma*j*delta).coefficient_tuple(), (gamma*j*i*delta).coefficient_tuple(), beta1.coefficient_tuple(), beta2.coefficient_tuple()]).transpose()
        C, D, _, _ = M.right_kernel().basis()[0]
    else:
        M = Matrix(Integers(I.norm()), [(gamma*j*delta).coefficient_tuple(), (gamma*j*i*delta).coefficient_tuple(), B(1).coefficient_tuple(), beta1.coefficient_tuple(), beta2.coefficient_tuple()]).transpose()
        C, D, _, _, _ = M.right_kernel().basis()[0]
    return ZZ(C), ZZ(D)

def FindLinearCombination(beta, I, K):
    r"""
    IdealToIsogenyEichler subroutine
    """
    i,j,k = I.quaternion_algebra().gens()
    alpha = IdealGenerator(I)
    twof = K.norm()
    O = I.left_order()
    assert I.left_order() == K.left_order()
    assert I + O*2 != K + O*2
    b1, b2, b3, b4 = K.gens() # Instead of using HNF, we can just add all the vectors...
    v_alpha = [c % twof for c in QuaternionOrderBasis(alpha, O)]
    v_abeta = [c % twof for c in QuaternionOrderBasis(alpha*beta, O)]
    v_b1 = [c % twof for c in QuaternionOrderBasis(alpha*b1, O)]
    v_b2 = [c % twof for c in QuaternionOrderBasis(alpha*b2, O)]
    v_b3 = [c % twof for c in QuaternionOrderBasis(alpha*b3, O)]
    v_b4 = [c % twof for c in QuaternionOrderBasis(alpha*b4, O)]

    def _kernel_mod_composite(A):
        m, n = A.dimensions()
        l = A.base_ring().characteristic()
        M = ZZ**m / (l*ZZ**m)
        phi = M.hom([(ZZ**n/(l*ZZ**n))(a) for a in A])
        for vec in phi.kernel().gens():
            vec_ker = M(vec)
            C, D, _, _, _, _ = vec_ker
            if gcd((C, D, 2)) == 1:
                return vec_ker
        assert False, "This should not happen if K is cyclic"

    M = Matrix(Integers(twof), [v_alpha, v_abeta, v_b1, v_b2, v_b3, v_b4])
    C, D, _, _, _, _ = _kernel_mod_composite(M)
    assert alpha*(C + D*beta) in K
    assert gcd([C, D, 2]) == 1
    return ZZ(C), ZZ(D)


#####################################
#                                   #
#         FullStrongApprox          #
#                                   #
#####################################

def FullStrongApproximation(O0, N_mu, N, C, D, lam, Kopt = None, alphaopt = None, q=1):
    r"""
    The "main" KLPT subroutine
    """
    B = O0.quaternion_algebra()
    i,j,k = B.gens()
    p = B.ramified_primes()[0]
    c0 = 2*p*lam*C % N
    c1 = pow(2*p*lam*D*q, -1, N)
    mu_0 = j*(C+i*D)
    N_mu0 = (mu_0).reduced_norm()
    c2 = ((N_mu - lam**2*N_mu0)/N) % N
    c3 = (-c0*c1) % N
    L = N*Matrix(ZZ, [[1, c3], [0, N]])
    v0 = vector(ZZ, [0, (c1*c2) % N])
    v = lam*vector(ZZ, [C, D]) + N*v0
    close_vectors = generate_close_vectors(L, -v, p, N_mu)
    for v_bar in close_vectors:
        z, w = v_bar + N*v0
        z, w = z/N, w/N
        M = ZZ(N_mu - p*((lam*C + z*N) + i*(lam*D + w*N)).reduced_norm())
        M = ZZ(M/(N**2))
        if M < 0:
            continue
        x, y, found = Cornacchia(O0.QF, M)
        if found:                
            mu_1 = x + i*y + j*z - k*w # -w because ji = -k
            mu = (lam*mu_0 + N*mu_1)
            if q == 1: #Only necessary in the secret cases
                x, y, z, w = mu.coefficient_tuple()
                if not (x % 2 == w % 2 and y % 2 == z % 2 and gcd([x,y,z,w]) == 1):
                    continue
            mu = MakePrimitive(mu, O0)
            if mu.reduced_norm() != N_mu//4:
                continue
            if alphaopt:
                mu = alphaopt*mu*alphaopt**(-1)
            passed = True 
            if Kopt:
                passed = (not mu in Kopt.left_order().intersection(Kopt.right_order()))
            if passed:
                return mu
    return None


#####################################
#                                   #
#           KLPT variants           #
#                                   #
#####################################

def KeyGenKLPT(O0, I, f):    
    r"""
    The KLPT version used in key generation
    """
    B = O0.quaternion_algebra()
    i,j,k = B.gens()
    p = B.ramified_primes()[0]
    N_I = ZZ(I.norm())
    assert N_I % 4 == 3
    Z_N = Integers(N_I)
    low_k = ceil(log(p,2) - log(N_I, 2) + 1)
    log_output_norm = ceil((3*log(p,2)+15)/f)*f #The +15 is guesstimated...
    for k in range(low_k, low_k+10):
        N_gamma = ZZ(2)**(k)
        N_mu = ZZ(2)**ZZ(log_output_norm - k)
        gamma = FullRepresentInteger(O0, N_I*N_gamma)
        if not gamma:
            continue
        C, D = EichlerModConstraint(I, gamma, 1, divisible=True)
        mu_0 = j*(C + i*D)
        lam_sqr = Z_N(4*N_mu)/Z_N(mu_0.reduced_norm())
        if not (lam_sqr).is_square():
            continue
        lam = ZZ(sqrt(lam_sqr))
        mu = FullStrongApproximation(O0, 4*N_mu, N_I, C, D, lam)
        if mu and gamma*mu == MakePrimitive(gamma*mu, O0):
            return gamma*mu
    return None

def SpecialEichlerNorm(T, O0, O0_alt, I, K, I_secret = None):
    r"""
    The KLPT version used in IdealToIsogenyEichler
    """
    B = O0.quaternion_algebra()
    p = B.ramified_primes()[0]
    if K:
        assert K.norm() == 2
        assert I.right_order() == K.left_order()
    if I_secret: #First step we must use the small secret ideal
        J = I_secret
        alpha = IdealEquivalence(I, J)
        beta = SpecialEichlerNormFixed(O0, T, J, K, alpha)
    else:
        if I.norm() < p**(0.4):
            J = None
        else:
            J = RandomEquivalentPrimeIdeal(I)
        if not J: #If the connecting ideal is unnaturally small, we need to use a different order
            B_new = O0_alt.quaternion_algebra()
            gamma, gamma_inv = IsomorphismGamma(B, B_new)
            O_target = B_new.quaternion_order([EvalIsomorphism(alpha, B_new, gamma) for alpha in I.right_order().gens()])
            I_new = ConnectingIdeal(O0_alt, O_target)
            if K:
                K_gen = IdealGenerator(K)
                K_new = O_target*EvalIsomorphism(K_gen, B_new, gamma) + O_target*K.norm()
            else:
                K_new = K
            O = O0_alt
        else: #If not, no need to be fancy.
            B_new = B
            gamma_inv = B(1)
            O = O0
            I_new = I
            K_new = K
        for _ in range(100):
            J = RandomEquivalentPrimeIdeal(I_new)
            if not J:
                continue
            alpha = IdealEquivalence(I_new, J)
            beta = SpecialEichlerNormFixed(O, T, J, K_new, alpha)
            if beta:
                beta = EvalIsomorphism(beta, B, gamma_inv)
                break
    if not beta:
        print("SpecialEichlerNorm failed miserably")
        assert False
    if K:
        assert not beta in K.left_order().intersection(K.right_order())
    assert beta in I.right_order()
    return beta

def SpecialEichlerNormFixed(O0, T, J, K, alpha):
    B = O0.quaternion_algebra()
    i,j,k = B.gens()
    p = B.ramified_primes()[0]
    q = abs(ZZ(i**2))

    N_J = ZZ(J.norm())
    Z_N = Integers(N_J)
    C, D = EichlerModConstraint(J, 1, 1)
    if C == 0 and D == 0:
        print(f"Output of EichlerModConstraint was {C, D}, which shouldnt happen")
        return None

    def _DivisorList(p, T, N_J):
        out = []
        for (l,e) in factor(T):
            if (T**2)/l > p*N_J**3:
                out.append(ZZ((T**2)/l))
        return out
    
    list_Nbeta = _DivisorList(p, T, N_J)
    assert len(list_Nbeta) > 0, "Ehmmm some parameters need to be increased"
    while len(list_Nbeta) > 0:
        N_beta = list_Nbeta.pop()
        mu_0 = j*(C + i*D)
        lam_sqr = Z_N(4*N_beta)/Z_N(mu_0.reduced_norm())
        if not (lam_sqr).is_square():
            continue
        lam = ZZ(sqrt(lam_sqr))
        beta = FullStrongApproximation(O0, 4*N_beta, N_J, C, D, lam, Kopt=K, alphaopt = alpha, q=q)
        if beta:
            return beta
    return None

def SigningKLPT(O0, K, I_tau, f):
    r"""
    The KLPT version used to compute the response ideal
    """
    B = O0.quaternion_algebra()
    i,j,k = B.gens()
    p = B.ramified_primes()[0]
    L = RandomEquivalentPrimeIdeal(K)
    if not L:
        return None
    Lm = pullback(I_tau, L)
    I = RandomEquivalentPrimeIdeal(Lm)
    if not I:
        return None
    delta = IdealEquivalence(Lm, I).conjugate()*Lm.norm()
    N_I = ZZ(I.norm())
    N_tau = ZZ(I_tau.norm())
    Z_N_I = Integers(N_I)
    Z_N_tau = Integers(N_tau)
    low_k = ceil(log(p,2) - log(N_I, 2) + 1)
    log_output_norm = ceil(((15/4)*log(p,2) + 25)/f)*f #Note for simplicity we make e an exact multiple of f
    for k in range(low_k, low_k+10):
        N_gamma = ZZ(2)**(k)
        N_mu = ZZ(2)**ZZ(log_output_norm - k)
        gamma = FullRepresentInteger(O0, N_I*N_gamma)
        if not gamma:
            continue
        C_0, D_0 = EichlerModConstraint(I, gamma, 1, divisible=True)
        a_0 = j*(C_0 + i*D_0)
        lam0_sqr = Z_N_I(N_mu)/Z_N_I(a_0.reduced_norm())
        if not (lam0_sqr).is_square():
            continue
        lam_0 = ZZ(sqrt(lam0_sqr))
        C_1, D_1 = EichlerModConstraint(I_tau, gamma, delta, divisible=False)
        a_1 = j*(C_1 + i*D_1)
        lam1_sqr = Z_N_tau(N_mu)/Z_N_tau(a_1.reduced_norm())
        if not (lam1_sqr).is_square():
            continue
        lam_1 = ZZ(sqrt(lam1_sqr))

        lam = CRT([lam_0, lam_1], [N_I, N_tau])
        C = CRT([C_0, C_1], [N_I, N_tau])
        D = CRT([D_0, D_1], [N_I, N_tau])
        mu = FullStrongApproximation(O0, 4*N_mu, N_I*N_tau, C, D, 2*lam)
        if not mu:
            continue
        beta = (gamma*mu*delta)/I.norm()
        if not (beta == MakePrimitive(beta, K.left_order())):
            continue
        if (beta).reduced_trace() % 2 == 1:
            return L*(beta.conjugate()/L.norm())
    return None
