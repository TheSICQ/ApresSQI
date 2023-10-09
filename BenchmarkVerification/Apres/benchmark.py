from ApresSQI import ApresSQI_verifier
from all_primes import all_primes
from fp_arith import update_p, reset_counter, get_cost
from tors_basis import basis_two_torsion, complete_basis_two_torsion_seed_get_plus, complete_basis_chall_torsion_seed_get_plus, point_difference
from ec_misc import is_supersingular, proj_point_normalize, double_proj_normalize, normalize_curve
from mont_xonly import ladder_3pt, xMULc, xMULaffs, xDBLe, xTPLe
from isogenychains import four_iso_chain_opt, three_iso_chain_strategy

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

def fake_sig_gen(verifier, A, variant, push = False):
    if variant in ['NIST', 'LWXZ', 'APRESnswp', 'APRESnsnp']:
        return gen_fake_no_seed(verifier)
    elif variant in ['APRESwsnp', 'APRESwswp']:
        return gen_fake_with_seed(verifier, A, push = push)
    elif variant == 'uncompressed':
        return gen_fake_uncompressed(verifier, A)
    else:
        assert False, "variant is wrong mate"

def gen_fake_no_seed(ver):
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

def gen_fake_with_seed(ver, A, push = False):
    zip, r, s = gen_fake_no_seed(ver)
    return get_seeds(ver, zip, r, s, A, push = push)

def gen_fake_uncompressed(ver, ProjA):
    zip, r, s = gen_fake_no_seed(ver)
    b, scalars = zip
    gens = []
    for count, si in enumerate(scalars):
        A = proj_point_normalize(ProjA)
        P, Q, _ = complete_basis_two_torsion_seed_get_plus(A, ver.f)
        P, Q = double_proj_normalize(P, Q)
        PmQ = point_difference(P, Q, A)        
        K = ladder_3pt(P, Q, PmQ, si, A)
        K = xMULaffs(K, ver.m, A)        
        if count != ver.B - 1 or ver.g == 0:
            ProjA, _ = four_iso_chain_opt(K, [], A, ver.f, ver.strategy_f)
        else:
            K = xDBLe(K, A, ver.f - ver.g)
            ProjA, _ = four_iso_chain_opt(K, [], A, ver.g, ver.strategy_g)
        gens.append(proj_point_normalize(K)[0])
    A = proj_point_normalize(ProjA)
    #(f"E_2 gen: {A}")
    P, Q, _ = complete_basis_two_torsion_seed_get_plus(A, ver.f2)
    if len(s) == 4:
        b, scalar, _, _ = s
    else:     
        b, scalar = s
    K = ladder_3pt(P, Q, PmQ, scalar, A)
    m2 = (ver.p+1)//(2**ver.f2)
    K = xMULaffs(K, m2, A)
    K, Q = double_proj_normalize(K, Q)
    ProjA, _ = four_iso_chain_opt(K, [], A, ver.f2, ver.strategy_f2)
    A, _ = normalize_curve(ProjA)
    #print(f"E_1 gen: {A}")
    if ver.f2 == 128:
        P, Q, hash_seeds = complete_basis_two_torsion_seed_get_plus(A, ver.f2)
    else:
        P, Q, hash_seeds = complete_basis_chall_torsion_seed_get_plus(A, ver.f2, ver.f3)
    sigma = (gens, A[0], hash_seeds, []) #forth empty list to indicate using uncompressed
    return sigma


def random_walk(ver, ProjA):
    #performs a small random walk in the 2-isogeny graph of a certain length and outputs A' in similar form after this walk
    P, Q, PmQ = basis_two_torsion(ProjA, ver.f2)
    s = randint(0, 2**128)
    K = ladder_3pt(P, Q, PmQ, s, ProjA)
    ProjA, Qlist = four_iso_chain_opt(K, [Q], ProjA, ver.f2, ver.strategy_f2)
    A = proj_point_normalize(ProjA)
    return A

def get_seeds(ver, zip, r, s, ProjA, push = False):
    starting_A = ProjA
    if push:
        seeds = [0 for _ in range(ver.B + 1)]
    else:
        seeds = [[0,0] for _ in range(ver.B + 1)]
    b, scalars = zip
    assert ver.B == len(scalars)
    for count, si in enumerate(scalars):
        A = proj_point_normalize(ProjA)
        if not push or count == 0:
            P, Q, seeds[count] = complete_basis_two_torsion_seed_get_plus(A, ver.f)
        else: 
            P, _, seed = complete_basis_two_torsion_seed_get_plus(A, ver.f)
            seeds[count] = seed[0]
        P, Q = double_proj_normalize(P, Q)
        #print(b)
        #print(f"gen {count}: {proj_point_normalize(A)}")
        #print(f"gen {count}: {proj_point_normalize(P), proj_point_normalize(Q)}")
        #print(f"Q above {count}: {proj_point_normalize(xMULaffs(Q, ver.m*2**(ver.f-1), A))}")
        PmQ = point_difference(P, Q, A)
        if b:
            temp = deepcopy(P)
            P = deepcopy(Q)
            Q = deepcopy(temp)
            b = 0
        K = ladder_3pt(P, Q, PmQ, si, A)
        K = xMULaffs(K, ver.m, A)
        #print(f"seeds: {seeds[count]}")
        #print(f"gen {count}: {proj_point_normalize(P), proj_point_normalize(Q)}")
        #print(f"K gen {count}: {proj_point_normalize(K)}")
        if count != ver.B - 1 or ver.g == 0:
            ProjA, Qlist = four_iso_chain_opt(K, [Q], A, ver.f, ver.strategy_f)
            Q = Qlist[0]
        else:
            K = xDBLe(K, A, ver.f - ver.g)
            ProjA, Qlist = four_iso_chain_opt(K, [Q], A, ver.g, ver.strategy_g)
            Q = Qlist[0]
    A = proj_point_normalize(ProjA)
    if ver.f3 == 0:
        if not push:
            P, Q, uchall = complete_basis_two_torsion_seed_get_plus(A, ver.f2)
        else:
            P, _, uchall = complete_basis_two_torsion_seed_get_plus(A, ver.f2)
            # Need to make sure its actually a basis
            uchall = uchall[0]
    else:
        if not push:
            P, Q, uchall = complete_basis_chall_torsion_seed_get_plus(A, ver.f2, ver.f3)
        else:
            P, _, uchall = complete_basis_chall_torsion_seed_get_plus(A, ver.f2, ver.f3)
            uchall = uchall[0]
            
    Q = proj_point_normalize(Q)
    PmQ = point_difference(P, Q, A)
    if len(s) == 2:
        b, scalar = s
        K = ladder_3pt(P, Q, PmQ, scalar, A)
        m2 = (ver.p+1)//(2**ver.f2)
        K = xMULaffs(K, m2, A)

        ProjA, _ = four_iso_chain_opt(K, [], A, ver.f2, ver.strategy_f2)

    # Combination of 2-, and 3-power
    elif len(s) == 4:    
        b1, s1, b2, s2 = s
        K2 = ladder_3pt(P, Q, PmQ, s1, A)
        if not b2:
            K3 = ladder_3pt(P, Q, PmQ, s2, A)
        else:
            K3 = ladder_3pt(Q, P, PmQ, s2, A)

        m2 = (ver.p+1)//(2**ver.f2)
        m3 = (ver.p+1)//(3**ver.f3)
        K2 = xMULaffs(K2, m2, A)
        K3 = xMULaffs(K3, m3, A)

        #Fix later: K3 might not have full order.
        if xTPLe(K3, A, ver.f3-1)[1] == [0, 0]: #Only happens in wswp
            return fake_sig_gen(ver, starting_A, 'APRESwswp', push=push)

        ProjA, Qlist = four_iso_chain_opt(K2, [K3], A, ver.f2, ver.strategy_f2)
        A = proj_point_normalize(ProjA)
        ProjA, Qlist = three_iso_chain_strategy(Qlist[0], [], A, ver.f3, ver.strategy_f3)
    ProjA, _ = normalize_curve(ProjA)
    #print(f"E_1 gen: {ProjA}")
    if ver.f3 == 0:
        P, Q, hash_seeds = complete_basis_two_torsion_seed_get_plus(ProjA, ver.f2)
    else:
        P, Q, hash_seeds = complete_basis_chall_torsion_seed_get_plus(ProjA, ver.f2, ver.f3)
    seeds_ls = [uchall, hash_seeds]
    return zip, r, s, seeds, seeds_ls

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
            variant = ApresVariants[i]
        plt.plot(x_coords, y_coords, label=variant)
        variant = None
    plt.ylim(0, 550000)
    plt.xlabel('f')
    plt.ylabel('cost')
    plt.legend()
    plt.show()

def cost_graph(variant, num_samples, primes, filename, doPlot = False):
    cost_array_samples = []

    ApresVariants = ['APRESwsnp', 'APRESwswp', 'APRESnswp', 'APRESnsnp', 'uncompressed']
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

        if f3 == 1 and variant == 'APRESwswp':
            continue #TODO fix this bug

        print(f'prime {p} with 2-torsion 2^{f}; D_chall = 2^{f2}*3^{f3}')
        assert((p+1) % (2**f2 * 3**f3) == 0)

        if variant in ApresVariants:
            SQI = ApresSQI_verifier(p, f, f2, f3)
            uncomp = False
            if variant == 'uncompressed':
                uncomp = True
            if 'wp' in variant:
                push = True
            else:
                push = False
        else:
            assert False, "variant is wrong mate"

        for _ in range(num_samples):
            #generate signature
            ProjA = random_walk(SQI, [[6,0],[1,0]])
            assert is_supersingular(ProjA)

            sigma = fake_sig_gen(SQI, ProjA, variant, push=push)

            reset_counter()
            #compute sig isogeny
            if variant in ApresVariants:
                _ = SQI.verify("Lets go for ApresSQI", sigma, ProjA, push = push)
            else:
                _ = SQI.verify("Lets go for ApresSQI", sigma, ProjA)
        
            cost_run = get_cost()
            #print([a - b for a, b in zip(counter_fp, counter_old)])
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
        if variant == 'uncompressed':
            primes = primes[1:]
    else:
        file = f'results_version_{variant}_samples_{num_samples}.txt'
    return cost_graph(variant, num_samples, primes, file, doPlot = plot)

ApresVariants = ['APRESwswp', 'APRESwsnp', 'APRESnswp', 'APRESnsnp', 'uncompressed']


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
        for variant in ApresVariants:
            data.append(benchmark_variant(variant, num_samples, plot = False, specific = specific))
        plot(data)
    else:
        benchmark_variant(variant, num_samples, plot = True, specific = specific)
