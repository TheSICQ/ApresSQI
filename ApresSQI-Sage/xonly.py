#!/usr/bin/env python3

from sage.all import *
from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite

class xPoint:
    r"""
    Class for x-only arithmetic on Montgomery curves.
    """
    def __init__(self, X, E):
        assert isMontgomery(E)
        k = E.base_field()
        self.X = k(X) if X is not None else None
        self.curve = E

    def __repr__(self):
        return f"Point with X-coordinate {self.X} on {self.curve}"

    def __bool__(self):
        return self.X is not None
    
    def __eq__(self, other):
        assert isinstance(other, xPoint)
        return self.X == other.X and self.curve == other.curve
    
    def xDBL(self):
        if not self:
            return self
        a = MontgomeryA(self.curve)
        X = self.X
        XX = X**2
        Z3 = 4*X*(XX + a*X + 1)
        if Z3 == 0:
            return xPoint(None, self.curve)
        X3 = (XX - 1)**2
        return xPoint(X3/Z3, self.curve)
    
    def push(self, isogeny):
        r"""
        Given an isogeny phi (where phi is computed as a composition of isogenies
        computed with Kohel's algorithm, or a composition of 2 and 4 isogenies), 
        returns phi(self)
        """
        if type(isogeny) is not EllipticCurveHom_composite:
            isogeny = EllipticCurveHom_composite.from_factors([isogeny])

        newX = self.X
        for isopart in isogeny.factors():
            assert isinstance(isopart, EllipticCurveIsogeny)
            if isopart._EllipticCurveIsogeny__algorithm == 'kohel':

                if (isom := isopart._EllipticCurveIsogeny__pre_isomorphism):
                    newX = isom.x_rational_map()(newX)

                phi = isopart._EllipticCurveIsogeny__phi
                psi = isopart._EllipticCurveIsogeny__psi
                try:
                    newX = phi(newX) / psi(newX)**2
                except ZeroDivisionError:
                    verbose("Point seems to be in the kernel")
                    newX = None
                    return xPoint(None, isogeny.codomain())

                if (isom := isopart._EllipticCurveIsogeny__post_isomorphism):
                    newX = isom.x_rational_map()(newX)

            elif 4 % isopart.degree() == 0:
                try:
                    newX = isopart.x_rational_map()(newX)
                except ZeroDivisionError:
                    verbose("Point seems to be in the kernel")
                    return xPoint(None, isogeny.codomain())
            else:
                print(isopart)
                print(isopart.degree())
                assert False, "Cannot evaluate in this type of isogeny"

        new_curve = isogeny.codomain().base_extend(self.curve.base_field())
        return xPoint(newX, new_curve)

    def xMUL(self, n):
        """
        Given an integer n, computes [n]self
        """
        n = ZZ(n).abs()
        if n == 0:
            return xPoint(None, self.curve)
        if n == 1:
            return self
        R0, R1, diff = self, self.xDBL(), self
        for i in [int(b) for b in bin(n)[3:]]:
            R0pR1 = xADD(R0, R1, diff)
            diff = xADD(R0, R1, R0pR1)
            if i == 0:
                R0, R1 = R0.xDBL(), R0pR1
            if i == 1:
                R0, R1 = R0pR1, R1.xDBL()
        return R0
    
    def x_multiples(self):
        r"""
        Given an xPoint P, 
        returns all the roots of the kernel polynomial of <P>.
        """
        if not self:
            return []
        Ps = [self]
        P = self.xDBL()
        if not P:
            xs = [P.X for P in Ps]
            return xs
        while P not in Ps[-2:]:
            Ps.append(P)
            P = xADD(Ps[-1], Ps[0], Ps[-2])
        xs = [P.X for P in Ps]
        return xs
    
    def kernel_polynomial(self, E, l):
        r"""
        Given that self has order l,
        returns the kernel polynomial of <P>
        """
        x = self.X
        F2 = E.base_field()
        Fbig = x.parent()
        ### wtf ###
        ext = Fbig.over(F2)
        ###########
        R,X = F2['X'].objgen()
        assert E.base_extend(Fbig) == self.curve

        assert l.is_prime()
        if l <= 3:
            return R([-x, 1])

        try:
            X_in_F2 = F2(x)
        except ValueError:
            pass
        else:
            return prod(X - xx for xx in xPoint(X_in_F2, E).x_multiples())

        if E.frobenius() not in ZZ:
            raise NotImplementedError

        def minpoly(elt):
            return ntl_funs(F2).minpoly_mod(ext(elt).vector(), ext.modulus())

        fs = [minpoly(x)]

        k = fs[0].degree()
        m = (l-1) // (2*k)

        assert k > 1    # handled above

        from sage.schemes.elliptic_curves.isogeny_small_degree import _least_semi_primitive
        a = _least_semi_primitive(l)

        xi = xPoint(x, E.change_ring(Fbig))
        for _ in range(1, m):
            xi = xi.xMUL(a)
            fs.append(minpoly(xi.X))

        return prod(fs)
    
    def xISOG(self, E, l):
        r"""
        Given a x = x(P) of a point P on E of order l,
        computes the separable isogeny with <P> as kernel
        """
        h = self.kernel_polynomial(E, l)
        phi = E.isogeny(h, check=False)
        _, phi = Normalized(phi)
        return phi

def xADD(P, Q, PmQ):
    r"""
    Given a x(P), x(Q) and x(P-Q),
    computes x(P+Q)
    """
    if not P:
        return Q
    if not Q:
        return P
    if not PmQ:
        return P.xDBL()
    A = P.X + 1
    B = P.X - 1
    C = Q.X + 1
    D = Q.X - 1
    DA = D*A
    CB = C*B
    X5 = (DA + CB)**2
    Z5 = PmQ.X*(DA - CB)**2
    if Z5 == 0:
        return xPoint(None, P.curve)
    return xPoint(X5/Z5, P.curve)

def xDBLMUL(m, n, P, Q, PmQ):
    r"""
    Algorithm 9 from eprint2017/212
    Given a x(P), x(Q) and x(P-Q), and scalars m,n
    computes x([m]P+[n]Q)
    """
    if m == 0:
        return Q.xMUL(n)
    if n == 0:
        return P.xMUL(m)
    s0, s1 = Integer(m), Integer(n)
    x0, x1, xmin = P, Q, PmQ
     
    while s0 != 0:

        if s1 < s0:
            s0, s1 = s1, s0
            x0, x1 = x1, x0

        if s1 <= 4*s0:
            s1 = s1 - s0
            tmp = xADD(x1, x0, xmin)
            xmin = x0
            x0 = tmp
        elif s0 % 2 == s1 % 2:
            s1 = (s1-s0)//2
            x0 = xADD(x1, x0, xmin)
            x1 = x1.xDBL()
        elif s1 % 2 == 0:
            s1 = s1//2
            xmin = xADD(x1, xmin, x0)
            x1 = x1.xDBL()
        else:
            s0 = s0//2
            xmin = xADD(x0, xmin, x1)
            x0 = x0.xDBL()

    while s1 % 2 == 0:
        s1 = s1//2
        x1 = x1.xDBL()

	#extra step with tripling as described in eprint2017/212
    while s1 % 3 == 0:
        s1 = s1//3
        x1 = x1.xMUL(3)

    if s1 > 1:
        x1 = x1.xMUL(s1)

    return x1

################################################################

import ast
class _ntl_funs:
    r"""
    An object encapsulating the NTL context for a given finite
    field F, such that polynomials in F[X] can be converted to
    NTL using .ntlify() and minimal polynomials in F[X]/f can
    be computed using .minpoly_mod().
    """
    def __init__(self, F):
        self.ctx = ntl.ZZ_pEContext(ntl.ZZ_pX(F.modulus().list(), F.characteristic()))
        self.F = F
        self.R = F['X']

    def ntlify(self, poly):
        try:
            poly = poly.vector()
        except AttributeError:
            pass
        assert poly.base_ring() == self.F
        return ntl.ZZ_pEX([ntl.ZZ_pE(c.list(), self.ctx) for c in poly])

    def minpoly_mod(self, elt, mod):
        ntl_mod = self.ntlify(mod)
        ntl_elt = self.ntlify(elt) % ntl_mod
        ntl_res = ntl_elt.minpoly_mod(ntl_mod)
        return self.R(ast.literal_eval(str(ntl_res).replace(' ',',')))

@cached_function
def ntl_funs(F2):
    r"""
    Caching helper function to construct the _ntl_funs object
    for a given base field F2 of degree 2 over its prime field.
    """
    assert F2.degree() == 2
    return _ntl_funs(F2)


#########################################
#
#   Unify 2 isogenies with verification
#
#########################################

def isMontgomery(E):
    r"""
    Given a curve E,
    returns whether E is in Montgomery form
    """
    a1, a2, a3, a4, a6 = E.a_invariants()
    return a1 == 0 and a3 == 0 and a6 == 0 and a4 == 1

def MontgomeryA(E):
    r"""
    Given a Montgomery curve E_A,
    returns the coefficient A
    """
    assert isMontgomery(E)
    return E.a2()

def projAC_toA(A24plus, C24):   
    A = 2*((2*A24plus) - C24)
    return A/C24

def two_iso_curve(P):
    r"""
    Given a point P of order 2, P != (0,0)
    returns the curve param A of the
    codomain of the isogeny generated by P
    """
    x = P.xy()[0]
    F = x.parent()
    assert x != 0
    A24plus = x**2
    C24 = F(1)
    A24plus = C24 - A24plus
    return projAC_toA(A24plus, C24)

def four_iso_curve(P):
    r"""
    Given a point P of order 4, [2]P != (0,0)
    returns the curve param A of the
    codomain of the isogeny generated by P
    """
    x = P.xy()[0]
    F = x.parent()
    C24 = F(1)
    A24plus = (1 + x**2)*(1 - x**2)
    return projAC_toA(A24plus, C24)

def four_iso_curve_singular(P, A):
    r"""
    Given a point P of order 4, [2]P = (0,0)
    returns the curve param A of the
    codomain of the isogeny generated by P
    """
    F = A.parent()
    K1 = F(4)
    ProjA24plus = (A + 2, K1)
    x = P.xy()[0]
    if x == 1:
        C24 = ProjA24plus[1] - ProjA24plus[0]
    elif x == -1:
        C24 = ProjA24plus[0]
    else:
        assert False, "something wrong with the point..."
    A24plus = K1
    return projAC_toA(A24plus, C24)

def four_iso(P, e):
    r"""
    Given a point P of order 2**e, [2**(e-1)]P != (0,0)
    returns the isogeny generated by [2**(e-2)]P
    """
    E_i = P.curve()
    F = E_i.base_field()
    P4 = 2**(e-2)*P
    next_A = four_iso_curve(P4)
    E_next = EllipticCurve(F, [0, next_A, 0, 1, 0])
    phi = E_i.isogeny(P4)
    psi = phi.codomain().isomorphism_to(E_next)
    return psi*phi

def customTwoIsogeny(P, f):
    r"""
    Given a point P of order 2**f,
    returns the isogeny generated by P
    """
    E_i = P.curve()
    F = E_i.base_field()
    assert isMontgomery(E_i)
    A = E_i.a2()

    # The first four isogeny
    P4 = 2**(f-2)*P
    if (2*P4).xy() == (0,0):
        next_A = four_iso_curve_singular(P4, A)
        phi = E_i.isogeny(P4)
    else:
        next_A = four_iso_curve(P4)
        phi = E_i.isogeny(P4)
    E_i = EllipticCurve(F, [0, next_A, 0, 1, 0])
    psi = phi.codomain().isomorphism_to(E_i)
    phi = psi * phi
    P = phi(P)
    f = f-2

    # The two iso, if f is odd
    if f%2 == 1:
        P2 = 2**(f-1)*P
        next_A = two_iso_curve(P2)
        phi_i = E_i.isogeny(P2)
        E_i = EllipticCurve(F, [0, next_A, 0, 1, 0])
        psi = phi_i.codomain().isomorphism_to(E_i)
        phi_i = psi*phi_i
        phi = phi_i * phi
        P = phi_i(P)
        f = f-1
    
    # The rest
    while f > 0:
        phi_i = four_iso(P, f)
        phi = phi_i*phi
        P = phi_i(P)
        f = f-2
    assert not P
    return phi

#############################
#                           #
#       Normalization       #
#                           #
#############################

def lexiographic_ordering(a):
    p = a.parent().characteristic()
    a0, a1 = list(a)
    return p*int(a0) + int(a1)

def sqrtDeterministic(a):
    r"""
    Returns a deterministic square root of a
    Very ugly to agree with fast square roots
    in signing, doing precomputation every time, but w/e
    """
    def get_constants(p):
        c0 = -1 # tries c0 = -1,-2,... for quadratic extension
        while pow(c0, (p-1)//2, p) == 1:
            c0 -= 1
        c0 = c0 % p  # smallest QNR in GF(p)
        c1 = pow(2, p-2 , p) # constant for division by 2  
        c2 = 1 # Tonelli-Shanks "e" in Michael Scott's paper
        while (p-1) % pow(2, c2 + 1, p) == 0:
            c2 += 1
        c3 = pow(c0, (p-1)*pow(pow(2, c2, p), p-2, p), p) # "z" in Michael Scott's paper
        c4 = (p-1-pow(2, c2, p))*pow(pow(2, c2 + 1, p), p-2, p) # y = x^c4 in Michael Scott's paper
        c4 = int(c4 % p)
        c6 = pow(2, c2, p) - 1
        c7 = c2 + 1
        c8 = pow(c0,p-2,p) #1/c0 in Michael Scott's paper
        c5 = pow(c0, c4, p)
        return [c0,c1,c2,c3,c4,c5,c6,c7,c8]
    
    def tonelli_shanks_update(a,b,cs):
        c0 = int(cs[0])
        c1 = cs[1]
        c2 = cs[2]
        c3 = cs[3]
        x = a
        y = b

        for _ in range(2, c0 + 1):
            x = x**2
            y = y**2

        y = y**2
        y = y*x

        if y != 1:
            a = a*c2
            b = b*c3

        x = a*b
        z = x*b

        for i in range(c0, 1, -1):
            w = z
            for _ in range(3, i+1):
                w = w**2
            if w != 1:
                x = x*c1
            c1 = c1**2
            if w != 1:
                z = z*c1
        return x, y
    
    F = a.parent()
    i = F.gens()[0]
    if i**2 != F(-1):
        return a.sqrt() #We only need to be consistent with signing
    p = F.characteristic()
    alist = list(a)
    c0,c1,c2,c3,c4,c5,c6,c7,c8 = get_constants(p)

    if a[1] == 0:
         if a[0].is_square():
              return F(a[0]**((p+1)/4))
         else:
              return i*a[0]**((p+1)/4)
         
    t0 = a[0]
    t1 = a[1]


    t2 = t0**2
    t3 = t1**2
    t3 = t3*c0
    t2 = t2-t3

    t3 = t2**c4

    t2, t3 = tonelli_shanks_update(t2, t3, [c2,c3,c0,c5])

    t3 = t0+t2
    t3 = t3*c1

    t0 = t3**c4

    t2 = t0

    for _ in range(1, c7 + 1):
        t2 = t2**2

    t2 = t2*c1
    t2 = t2*t1

    t1 = t3**c6

    t2 = t1*t2

    t3, t0 = tonelli_shanks_update(t3, t0, [c2,c3,c0,c5]) 

    t2 = t2*t3

    if t0 == 1:
        t0 = t3
        t3 = t2
    else:
        t0 = t2
        t3 = t3*c8

    a_sqrt = F([t0, t3])



    if int(t0) % 2 != 0: ## TODO: understand why this turns out this is the correct choice of squareroot
        return F([t0*c0, t3*c0]) 
    return F([t0, t3])

def MontgomeryNormalize(A):
    r"""
    Given a Montgomery coefficient A,
    returns a deterministic choice of
    Montgomery coefficient for the 
    isomorphism class [E_A]
    """
    F = A.parent()
    Z_0 = A**2
    temp = (A**3 - 3*A)/(2*sqrtDeterministic(A**2 - 4))
    Z_1 = (9 - A**2)/2 + temp 
    Z_2 = (9 - A**2)/2 - temp
    Zlist = [Z_0, Z_1, Z_2]
    Zlist.sort(key = lexiographic_ordering)
    Z = Zlist[0]
    Anew = sqrtDeterministic(Z)
    E = EllipticCurve(F, [0,Anew,0,1,0])
    return E

def Normalized(phi, R=None):
    r"""
    Given an isogeny phi,
    returns phi composed with an isomorphism
    to a deterministic choice of
    Montgomery coefficient for the 
    isomorphism class of the codomain.
    Additionally, given R, it also returns
    the evaluation of R.
    """
    A = phi.codomain().montgomery_model().a2()
    E_A = MontgomeryNormalize(A)
    psi = phi.codomain().isomorphism_to(E_A)
    if R:
        return E_A, psi*phi, psi(R)
    return E_A, psi*phi
        
