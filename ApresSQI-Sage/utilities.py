
####################################################################################################################
####################################  This file contains miscellaneous functions  ##################################
####################################################################################################################

#####################################################
#                                                   #
#                  from LearningToSQI               #
# https://github.com/LearningToSQI/SQISign-SageMath #
#                                                   #
#####################################################

"""MIT License

Copyright (c) 2023 Maria Corte-Real Santos, Jonathan Komada Eriksen, Michael Meyer  and Giacomo Pope

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software."""

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