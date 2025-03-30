
####################################################################################################################
#####################################         This file is used for testing      ###################################
####################################################################################################################


from klpt import *
from id2iso import *
from ApresSQI import *
from xonly import *
from xonly_velusqrt import EllipticCurveHom_velusqrt
import time

def test_RepresentInteger(param):
    print("Testing RepresentInteger")
    params = SQIsign(param)
    M = params.p*123 + 1
    O0 = params.O0
    gamma = FullRepresentInteger(O0, M)
    assert gamma in O0
    assert not gamma/2 in O0
    assert M == gamma.reduced_norm()
    print("  > Success")

def test_SimpleTranslations(param):
    print("Testing simple ideal-isogeny translations of a T-isogeny")
    params = SQIsign(param)
    alpha_test = FullRepresentInteger(params.O0, params.T*next_prime(params.p))
    I_test = params.O0*alpha_test + params.O0*params.T
    t_start = time.time()
    phi_gens = IdealToIsogenyGens(params.O0, I_test, params.facToBasis, params.facToAction)
    phi = chain_iso(phi_gens, params.E0)
    print(f"Translation finished: Took a total of {time.time() - t_start} seconds")
    I_testnew = IsogenyGensToIdeal(params.O0, phi_gens, params.facToAction)
    assert I_testnew == I_test
    print("  > Success")

def test_KeyGenKLPT(param):
    print("Testing KeyGenKLPT")
    params = SQIsign(param)
    print("Precomputation done")        
    p_bitsize = ceil(log(params.p,2))
    I_secret_bitsize = ceil(p_bitsize/4) 
    # The secret degree
    D_secret = next_prime(randint(2**I_secret_bitsize, 2**(I_secret_bitsize+1)))
    while D_secret % 4 != 3:
        D_secret = next_prime(D_secret)
    gamma = FullRepresentInteger(params.O0, D_secret * 2**p_bitsize)
    i,_,_ = params.B.gens()
    # The alternative path
    print("Starting KLPT")
    for tries in range(100):
        print(f"try: {tries}")
        a = randint(0, D_secret)
        I_secret = params.O0*(gamma*(a + i)) + params.O0*D_secret
        alpha = KeyGenKLPT(params.O0, I_secret, params.f)
        if alpha:
            break
    assert alpha in I_secret
    J_secret = I_secret*(alpha.conjugate()/D_secret)
    assert all([l == 2 for l,e in factor(J_secret.norm())])
    print("  > Success")

def test_SigningKLPT(param):
    print("Testing SigningKLPT")  
    params = SQIsign(param)
    print("Precomputation done")

    p_bitsize = ceil(log(params.p,2))
    I_secret_bitsize = ceil(p_bitsize/4) 
    # The secret degree
    D_secret = next_prime(randint(2**I_secret_bitsize, 2**(I_secret_bitsize+1)))
    while D_secret % 4 != 3:
        D_secret = next_prime(D_secret)
    # The secret ideal
    tau = FullRepresentInteger(params.O0, D_secret * 2**p_bitsize)
    I_tau = params.O0*tau + params.O0*D_secret
    cc = FullRepresentInteger(params.O0, params.T*(2**params.f)*next_prime(D_secret))
    I_comchall = params.O0*cc + params.O0*params.T**2*(2**params.f)
    K = I_tau.conjugate()*I_comchall
    for _ in range(100):
        J = SigningKLPT(params.O0, K, I_tau, params.f)
        if J:
            break
    assert IdealEquivalence(J, K)
    assert J.norm().radical() == 2
    print("  > Success")


def test_SpecialEichlerNorm(param):
    print("Testing SpecialEichlerNorm")
    params = SQIsign(param)
    print("Precomputation done")
    tau = FullRepresentInteger(params.O0, 2**params.f*next_prime(params.p))
    I_tau = params.O0*tau + params.O0*2**params.f
    beta = SpecialEichlerNorm(params.T, params.O0, params.O0_alt, I_tau, None)
    assert beta in I_tau.right_order()
    assert params.T**2 % beta.reduced_norm() == 0
    print("  > Success")

def test_IdealToIsogenyEichler(param):
    print("Testing IdealToIsogenyEichler")
    params = SQIsign(param)
    g = params.f*5
    tau = FullRepresentInteger(params.O0, 2**g*params.T)
    I_tau = params.O0*tau + params.O0*(2**g)
    phi_Itau, zip, _ = IdealToIsogenyEichler(params.O0, params.O0_alt, I_tau, params.O0*1, params.facToBasis, params.facToAction, params.B_2[0], params.f, params.T)
    I_alt = params.O0*tau.conjugate() + params.O0*(params.T)
    phi_alt = IdealToIsogeny(params.O0, I_alt, params.E0, params.facToBasis, params.facToAction)
    assert phi_Itau.codomain().j_invariant() == phi_alt.codomain().j_invariant()
    print("  > Success")

def test_KeyGen(param):
    print("Testing Key generation (KeyGen)")
    t_start = time.time()
    signer = SQIsign(param)
    pk = signer.KeyGen()
    print(f"Key generation finished: Took a total of {time.time() - t_start} seconds")
    signer.sk.export('SQI_' + param + '_priv.key')
    print("  > Success")

def test_Signing(param, compressed = True):
    print("Testing Signing")
    signer = SQIsign(param)
    #signer.KeyGen()
    signer.load_privkey(param=param)
    sigma = signer.Sign("Øl, øl, øl og aftersqi!", compressed = compressed)
    print(sigma)
    print("  > Success")

def test_SigningAndVerif(param, seeded = True, compressed = True):
    print("Testing Signing and verification together")
    signer = SQIsign(param)
    #signer.KeyGen()
    signer.load_privkey(param=param)
    msg = "Øl, øl, øl og aftersqi!"
    t_start = time.time()
    sigma = signer.Sign(msg, seeded = seeded, compressed = compressed)
    print(f"Signing finished: Took a total of {time.time() - t_start} seconds")
    print(sigma)
    ("Verifying...")
    verified = signer.verify(msg, sigma, signer.pk)
    print('Signature was ' + ('CORRECT' if verified else 'WRONG'))
    assert verified
    print("  > Success")

def test_xonly():
    print("Testing xOnly arithmetic")
    p = 507227047723007
    F2 = GF((p,2), name='z2', modulus=var('x')**2 + 1)   
    A = F2([62261760287349,159667122451046])
    E = EllipticCurve(F2, [0, A, 0, 1, 0])
    P, Q = E.gens()
    PmQ = P - Q
    a, b = [randint(0,1000000) for _ in range(2)]
    xP = xPoint(P.xy()[0], E)
    xQ = xPoint(Q.xy()[0], E)
    xPmQ = xPoint(PmQ.xy()[0], E)
    control = a*P + b*Q
    xPmQ = xPoint(PmQ.xy()[0], E)
    res = xDBLMUL(a, b, xP, xQ, xPmQ)
    assert res.X == control.xy()[0]
    print("  > Success")

def test_xSqrtVelu():
    p = 13632396072050491391
    F2 = GF((p,2), name='z2', modulus=var('x')**2 + 1)   
    E = EllipticCurve(F2, [1, 0])
    K, _ = E.gens()
    l = 1009
    K *= sqrt(E.order())//l
    E, _ = Normalized(E.isogeny(K))

    K, _ = E.gens()
    l = 1009
    K *= sqrt(E.order())//l

    xK = xPoint(K.xy()[0], E)
    phi = EllipticCurveHom_velusqrt(E, xK, l)

    phi_control = E.isogeny(K)
    phi_control = phi_control.codomain().isomorphism_to(phi.codomain()) * phi_control

    P = E.random_point()
    print(phi_control(P))
    print(phi(P))
    print(xPoint(P.xy()[0], E).push(phi))

    print("  > Success")

def CGL_collision(param):    
    params = SQIsign(param)    
    
    g = params.f

    while g < 1024:
        g += params.f

    #gamma = FullRepresentInteger(params.O0, Integer(2)**(2*g))
    i, j, k = params.B.gens()
    gamma = -Integer(716699723090500146378365648476854567893141048600647231835968645114838313947539777283490639034351049315596689297466352028430557117882871659176647956068067651829181425340409180954406537963130503832225155132673823429361752561773661442340566883844329563438846929162058132940359899435753233377244054329264015876212204635201) + Integer(574404203146155031533094788205053116156580814071514608434543574548868243857511193455582660605175872467306345807348830847196027495297476129881395716644629191082940341344567411930012967182006312446898356713459483142574462106046407746776304408435032505742907387768224625329922927453621307044856820895301835394305133789095)/2*i - Integer(7)/2*j - Integer(7)*k
    print(gamma)
    J1 = params.O0*gamma + params.O0*Integer(2)**(g)
    J2 = params.O0*gamma.conjugate() + params.O0*Integer(2)**(g)
    Js = [J1, J2]

    I_secret = RandomEquivalentPrimeIdeal(J1)
    while I_secret.norm() % 4 != 3:
        I_secret = RandomEquivalentPrimeIdeal(J1)
    while len(Js) < 3:
        print("Starting KLPT")
        for tries in range(100):
            print(f"try: {tries}")
            alpha = KeyGenKLPT(params.O0, I_secret, params.f)
            if alpha:
                break
        assert alpha in I_secret
        J_secret = I_secret*(alpha.conjugate()/I_secret.norm())

        assert all([l == 2 for l,e in factor(J_secret.norm())])
        Js.append(J_secret)

    for J_id in Js:
        phi_Itau, zip, _ = IdealToIsogenyEichler(params.O0, params.O0_alt, J_id, params.O0*1, params.facToBasis, params.facToAction, params.B_2[0], params.f, params.T)
        print("!!!!! PATH:")
        for phi in phi_Itau.factors():
            print(phi.degree(), phi.codomain().j_invariant())
        print("  > Success")

def all_tests(param):
    test_RepresentInteger(param)
    test_KeyGenKLPT(param)
    test_SpecialEichlerNorm(param)
    test_SigningKLPT(param)
    test_SimpleTranslations(param)
    test_IdealToIsogenyEichler(param)
    test_SigningAndVerif(param)
    test_SigningAndVerif(param, compressed=False)
    test_SigningAndVerif(param, seeded=False)

if __name__=="__main__":
    #param = '4-block'
    #param = '7-block'
    #param = 'NIST'
    #param = 'toy'
    #param = '136319'
    #param = '8513034219037441780170691209753296498696014329521974009944792576819199999999'
    param = '73743043621499797449074820543863456997944695372324032511999999999999999999999'
    #test_KeyGen(param)
    #test_SigningAndVerif(param)
    #test_SigningAndVerif(param, compressed=False)
    #test_SigningAndVerif(param, seeded=False)
    #all_tests(param)
    #test_xSqrtVelu()
    #test_IdealToIsogenyEichler(param)
    #test_KeyGen(param)
    #test_Signing(param)
    CGL_collision(param)
    print('All tests passed!')
