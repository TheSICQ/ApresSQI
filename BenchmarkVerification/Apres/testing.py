from ApresSQI import *

def valuation(n, ell):
    res = 0
    while not n % ell:
        res += 1
        n //= ell
    return res

def test_verif_uncompressed():
    p = 507227047723007
    update_p(p)
    sig = ([[444523460281294, 78583382474725], [132898737859456, 384888631179922], [322716351208294, 126891443424223], [187893610113812, 185273325662937], [124652030358258, 142334843532270], [421219858565922, 49097314478060], [350980029151906, 385739757700056]], [212335589880610, 142648608012200], (3, 5), [])
    f = valuation(p+1, 2)
    f3 = 0
    verifier = ApresSQI_verifier(p, f, f, f3)
    #This real example uses full size blocks
    verifier.e = ceil(((15/4)*log(verifier.p,2) + 25)/verifier.f)*verifier.f
    verifier.B = ceil(verifier.e/verifier.f)
    verifier.g = 0
    #40861744919524*z2 + 108841877744184
    projA = [[108841877744184,40861744919524],[1,0]]
    msg = "Øl, øl, øl og aftersqi!"
    verified = verifier.verify(msg, sig, projA)
    print(f'Uncompressed signature on msg = "{msg}" was ' + ('CORRECT' if verified else 'WRONG'))
    assert verified
    print("    > DONE")

def test_verif_wsnp():
    p = 507227047723007
    update_p(p)
    sig = ((1, [5014353738, 650328519, 3966955375, 678881623, 2438798155, 302252207, 3031629614]), 2541095627, (0, 2528748734), [(0, 3), (1, 1), (1, 0), (0, 5), (0, 1), (0, 1), (0, 0)], [(0, 1), (1, 7)])
    f = valuation(p+1, 2)
    f3 = 0
    verifier = ApresSQI_verifier(p, f, f, f3)
    #This real example uses full size blocks
    verifier.e = ceil(((15/4)*log(verifier.p,2) + 25)/verifier.f)*verifier.f
    verifier.B = ceil(verifier.e/verifier.f)
    verifier.g = 0
    #40861744919524*z2 + 108841877744184
    projA = [[108841877744184,40861744919524],[1,0]]
    msg = "Øl, øl, øl og aftersqi!"
    verified = verifier.verify(msg, sig, projA, push=False)
    print(f'Seeded signature on msg = "{msg}" was ' + ('CORRECT' if verified else 'WRONG'))
    assert verified
    print("    > DONE")


if __name__ == "__main__":
    #test_verification()
    #test_verification23()
    test_verif_uncompressed()
    test_verif_wsnp()
