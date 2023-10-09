from SQIsign import Verifier
from all_primes import all_primes
from fp_arith import update_p, reset_counter, get_cost
from tors_basis import basis_two_torsion
from ec_misc import is_supersingular, proj_point_normalize
from mont_xonly import ladder_3pt
from isogenychains import four_iso_chain_nist

from copy import deepcopy
import matplotlib.pyplot as plt
from math import ceil
from random import randint
import sys

def valuation(n, p):
    res = 0
    while not n % p:
        res += 1
        n //= p
    return res

def fake_sig_gen(ver, A, variant):
    assert variant in ['NIST', 'LWXZ']
    assert 2**ver.f2*3**ver.f3 >= 2**128
    #b = randint(0,1) TODO: Fix bug when this is 1...
    b = 0
    slist = []
    for _ in range(ceil(ver.e/ver.f)):
        slist.append(randint(1, 2**ver.f-1))
    zip = (b, slist)
    r = randint(0, 2**128)
    s = [0, randint(1, ver.f2)] 
    if ver.f3 > 0:
        s.append(0) #TODO: Fix bug when this is 1... (though it shouldnt matter for benchmarking)
        s.append(randint(1, 3**ver.f3))
    return (zip, r, s)

def random_walk(ver, ProjA):
    #performs a small random walk in the 2-isogeny graph of a certain length and outputs A' in similar form after this walk
    P, Q, PmQ = basis_two_torsion(ProjA, ver.f2)
    s = randint(0, 2**128)
    K = ladder_3pt(P, Q, PmQ, s, ProjA)
    ProjA, _ = four_iso_chain_nist(K, [], ProjA, ver.f2, ver.strategy_f2)
    A = proj_point_normalize(ProjA)
    return A 

def get_averages(cost_array):
    output = {}
    for f, cost in cost_array:
        if f in output.keys():
            output[f].append(cost)
        else:
            output[f] = [cost]
    return [(f, sum(costlist)/len(costlist)) for f, costlist in output.items()]

def plot(data, variant = None):
    for i, coords in enumerate(data):
        x_coords, y_coords = zip(*coords)
        if not variant:
            variant = Variants[i]
        plt.plot(x_coords, y_coords, label=variant)
        variant = None
    plt.ylim(0, 550000)
    plt.xlabel('f')
    plt.ylabel('cost')
    plt.legend()
    plt.show()

def cost_graph(variant, num_samples, primes, filename, doPlot = False):
    cost_array_samples = []

    for i, p in enumerate(primes):
        update_p(p)
        #set params
        f = valuation(p+1, 2)
        if f < 128:
            f2 = f
            f3 = valuation(p+1, 3)
        else:
            f2 = 128
            f3 = 0

        print(f'prime {p} with 2-torsion 2^{f}; D_chall = 2^{f2}*3^{f3}')
        assert((p+1) % (2**f2 * 3**f3) == 0)

        SQI = Verifier(p, f, f2, f3, version = variant)
        for _ in range(num_samples):
            #generate signature
            ProjA = random_walk(SQI, [[6,0],[1,0]])
            assert is_supersingular(ProjA)

            sigma = fake_sig_gen(SQI, ProjA, variant)

            reset_counter()
            #compute sig isogeny
            _ = SQI.verify("Lets go for ApresSQI", sigma, ProjA)
        
            cost_run = get_cost()
            print(f'total cost: {cost_run}')
            cost_array_samples.append([f, cost_run])
    cost_array = get_averages(cost_array_samples)
    if filename:
        with open(filename, "w") as file:
            file.write(str(cost_array))
    if doPlot:
        plot([cost_array])
    return cost_array

def benchmark_variant(v, num_samples, plot = False, specific = False):
    primes = all_primes
    if specific:
        primes = SpecificPrimes
        file = f'results_specific_version_{variant}.txt'
    else:
        file = f'results_version_{variant}_samples_{num_samples}.txt'
    return cost_graph(variant, num_samples, primes, file, doPlot = plot)

Variants = ['NIST', 'LWXZ']

SpecificPrimes = [23920667128620486487914848107166358953830561597426178123910317653495243603967, 21986677567972288250995905822739208616445482086236719868210978703712341458943, 22728720641309136015759539049556903787604752849407962277276342173428260798463]

if __name__ == "__main__":
    variant = 'All'
    if len(sys.argv) > 1:
        variant = sys.argv[1]

    num_samples = 5
    if len(sys.argv) > 2:
        num_samples = int(sys.argv[2])
    
    specific = False
    if len(sys.argv) > 3:
        specific = True
        print("running for specific primes!")

    logging = True

    #if logging:
        #print(f'logging is turned ON, results are written to file results_version_{variant}_samples_{num_samples}.txt\n\n')
    #else:
        #print("logging is turned OFF, results are not written to file\n\n")
    if variant == 'All':
        data = []
        for variant in Variants:
            data.append(benchmark_variant(variant, num_samples, plot = False, specific = specific))
        plot(data)
    else:
        benchmark_variant(variant, num_samples, plot = True, specific = specific)
