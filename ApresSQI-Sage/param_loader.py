from sage.all import *
from ast import literal_eval
from xonly import xPoint

def load(p):
    if p == 'toy':
        p = 507227047723007
    if p == 'toy_even':
        p = 117647744172031
    if p == 'NIST':
        p = 0x34e29e286b95d98c33a6a86587407437252c9e49355147ffffffffffffffffff
    if p == '4-block':
        p = 0x323fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
    if p == '7-block':
        p = 0x309c04bcaedbb0134cca8373e439ffffffffffffffffffffffffffffffffffff
    params = {}
    filename = 'Precomputed/' + str(p) + '.py'
    print(filename)
    try:
        with open(filename, "r") as file:
            p = Integer(literal_eval(file.readline()))
            f = Integer(literal_eval(file.readline()))
            D_com = Integer(literal_eval(file.readline()))
            D_chall = Integer(literal_eval(file.readline()))
            T = Integer(literal_eval(file.readline()))
            B_2 = literal_eval(file.readline())
            B_chall = literal_eval(file.readline())
            facToExt = literal_eval(file.readline())
            degToMod = literal_eval(file.readline())
            facToBasis_str = literal_eval(file.readline())
            facToAction_list = literal_eval(file.readline())
    except:
        print(f"file {filename} not found!")
        assert False, 'not precomputed yet!'
    
    B = QuaternionAlgebra(-1, -p)
    i, j, k = B.gens()
    O0 = B.quaternion_order([1, i, (i+j)/2, (1+k)/2])

    F = GF((p,2), name='z2', modulus=var('x')**2 + 1)    
    E0 = EllipticCurve(F, [1,0])
    E0.set_order((p+1)**2, num_checks=0)
    sqrtm1 = F(-1).sqrt()

    degToField = {1 : F}
    Zx, x = PolynomialRing(ZZ, 'x').objgen()
    for extdeg in degToMod.keys():
        Fbig = GF((p,2*extdeg), name='z' + str(2*extdeg), modulus=Zx(degToMod[extdeg]))
        if not Fbig.has_coerce_map_from(F):
            Fbig.register_coercion(Hom(F, Fbig).an_element())
        degToField[extdeg] = Fbig

    B_2 = [E0([F(c) for c in P]) for P in B_2]
    B_chall = [E0([F(c) for c in P]) for P in B_chall]

    facToBasis = {}
    for fac in facToBasis_str.keys():
        if fac%2==0: #either 2 torsion or D_chall torsion
            extdeg = 1
            Fbig = F
        else:
            extdeg = facToExt[fac]
            Fbig = degToField[extdeg]
        Ebig = E0.base_extend(Fbig)
        order = p**extdeg - (-1)**extdeg
        Ebig.set_order(order**2, num_checks=0)
        P, Q, PmQ = [xPoint(Fbig(P), Ebig) for P in facToBasis_str[fac]]
        facToBasis[fac] = [P, Q, PmQ]

    P2, Q2 = B_2
    xB_2 = xPoint(P2.xy()[0], E0), xPoint(Q2.xy()[0], E0), xPoint((P2-Q2).xy()[0], E0)
    facToBasis[2**f] = xB_2

    Pchall, Qchall = B_chall
    xB_chall = xPoint(Pchall.xy()[0], E0), xPoint(Qchall.xy()[0], E0), xPoint((Pchall-Qchall).xy()[0], E0)
    facToBasis[D_chall] = xB_chall

    facToAction = {}
    for fac in facToAction_list.keys():
        facToAction[fac] = [Matrix(Integers(fac), M) for M in facToAction_list[fac]]

    ### ApresSQI special: Precomputed values for use in basis generation
    z2 = F.gens()[0]
    small_nonsquares = []
    for n in range(10000):
        if not (n + z2).is_square():
            small_nonsquares.append(n+z2)
        if len(small_nonsquares) > 2**8:
            break
    small_square_min_one_nonsquare = []
    for n in range(1, 10000):
        if (n + z2).is_square() and not (n - 1 + z2).is_square():
            small_square_min_one_nonsquare.append(n+z2)
        if len(small_square_min_one_nonsquare) > 2**8:
            break

    params['p'] = p
    params['F'] = F
    params['f'] = f
    params['D_com'] = D_com
    params['D_chall'] = D_chall
    params['T'] = T
    params['B_2'] = B_2
    params['B_chall'] = B_chall
    params['facToExt'] = facToExt
    params['facToBasis'] = facToBasis
    params['facToAction'] = facToAction
    params['small_ns'] = small_nonsquares
    params['small_s'] = small_square_min_one_nonsquare

    return params

def AlternativeOrder(p):
    r"""
    Generates another special p-extremal maximal order.
    Sometimes useful in SpecialEichlerNorm
    """
    def _alt_QAlg(p):
        q = 1
        while True:
            q = next_prime(q)
            if (-q) % 4 == 1 and kronecker(-q, p) == -1:
                return QuaternionAlgebra(-q, -p), q
        
    B_alt, q = _alt_QAlg(p)
    O0_alt = B_alt.maximal_order()
    O0_alt.QF = BinaryQF([1,0,q])
    return O0_alt