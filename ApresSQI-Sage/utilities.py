from sage.all import *
import cypari2
pari = cypari2.Pari()

def discrete_log_pari(a, base, order):
    r"""
    Wrapper around pari discrete log.
    Works like a.log(b), but allows
    us to use the optional argument
    order.
    """
    x = pari.fflog(a, base, order)
    return ZZ(x)

def BiDLP(R, P, Q, D): 
    r"""
    Given points R, and a D-torsion basis P, Q
    returns a, b s.t. R = [a]P + [b]Q
    """
    ePQ = P.weil_pairing(Q, D, algorithm="pari")
    eRQ = R.weil_pairing(Q, D, algorithm="pari")
    eRP = R.weil_pairing(-P, D, algorithm="pari")

    a = discrete_log_pari(eRQ, ePQ, D)
    b = discrete_log_pari(eRP, ePQ, D)
    assert R == a*P + b*Q
    return a, b