from strategies_sike import *
from isogenychains import *
from copy import deepcopy
import matplotlib.pyplot as plt
import sys
from math import ceil
from NIST_SQIsign import NIST_verifier
from LWXZ_SQIsign import LWXZ_verifier
from APRES_wsnp_SQIsign import APRES_wsnp_verifier
from APRES_wswp_SQIsign import APRES_wswp_verifier
from APRES_nswp_SQIsign import APRES_nswp_verifier
from APRES_nsnp_SQIsign import APRES_nsnp_verifier

def valuation(n, p):
    res = 0
    while not n % p:
        res += 1
        n //= p
    return res

def fake_sig_gen(f2, f3):
        assert 2**f2*3**f3 >= 2**128
        b = randint(0,1)
        slist = []
        for _ in range(ceil(e/f2)):
            slist.append(randint(1, 2**f2-1))
        zip = (b, slist)
        r = randint(0, 2**128)
        s = [randint(0, 1), randint(1, 2**min(128, f2))]
        if f3 > 0:
            s.append(randint(0, 1))
            s.append(randint(1, 3**f3 - 1))
        return (zip, r, s)     

def graph_from_data(data):

    x_coords, y_coords = zip(*data)
    plt.plot(x_coords, y_coords, label='Connected Line', color='red')
    plt.xlabel('f')
    plt.ylabel('cost')
    
    return plt

def sum_f2(cost_array):
    tmplist = list(set(id for id, tmp in cost_array))
    tmplist.sort()
    return [
    [id, sum(val for _, val in cost_array if _ == id) // len([val for _, val in cost_array if _ == id])]
    for id in tmplist
    ]


def cost_graph(variant, num_samples, primes, file):
    cost_array_samples = []

    #go through all primes defined in globals.py
    for i, p in enumerate(primes):
        # if i%2 == 0:
        #     continue

        update_p(p)

        #set params
        #f2 = Integer(p + 1).valuation(2)
        f2 = valuation(p+1, 2)
        if f2 < 128:
            #f3 = Integer(p + 1).valuation(3)
            f3 = valuation(p+1, 3)
        else:
            f3 = 0

        print(f'prime p[{i}] with 2-torsion {f2} and 3-torsion {f3}')
        assert((p+1) % (2**f2 * 3**f3) == 0)

        if variant == 'NIST':
            SQI = NIST_verifier(p, f2, f3)
        elif variant == 'LWXZ':
            SQI = LWXZ_verifier(p, f2, f3)
        elif variant == 'APRESwsnp':
            SQI = APRES_wsnp_verifier(p, f2, f3)        
        elif variant == 'APRESwswp':
            SQI = APRES_wswp_verifier(p, f2, f3)        
        elif variant == 'APRESnswp':
            SQI = APRES_nswp_verifier(p, f2, f3)
        elif variant == 'APRESnsnp':
            SQI = APRES_nsnp_verifier(p, f2, f3)
        else:
            print("variant is wrong mate")
            return 0
        
        for _ in range(num_samples):
            #generate signature
            sigma = fake_sig_gen(f2, f3)
            ProjA = SQI.random_walk([[6,0],[1,0]])

            assert is_supersingular(ProjA)

            if variant == 'APRESwsnp' or variant == 'APRESwswp':
                ProjB = ProjA.copy()
                sigma = SQI.compute_seeds("Lets go for ApresSQI", sigma, ProjB)

            counter_old = deepcopy(counter_fp)
            #compute sig isogeny
            _ = SQI.verify("Lets go for ApresSQI", sigma, ProjA)
        
            cost_run = cost([a - b for a, b in zip(counter_fp, counter_old)], metric_fp)
            #print([a - b for a, b in zip(counter_fp, counter_old)])
            print(f'total cost: {cost_run}')
            cost_array_samples.append([f2, cost_run])
        
        if logging:
            append_to_file(file, f"{f2}, {sum_f2(cost_array_samples)[-1][1]}\n")


    cost_array = sum_f2(cost_array_samples)

    # plot = graph_from_data(cost_array)
    # plot.show()

    return cost_array_samples

def append_to_file(filename, content):
    with open(filename, 'a') as file:
        file.write(content)

if __name__ == "__main__":

    variant = 'NIST'
    if len(sys.argv) > 1:
        variant = sys.argv[1]

    num_samples = 1
    if len(sys.argv) > 2:
        num_samples = int(sys.argv[2])

    start_primes = 0
    if len(sys.argv) > 3:
        start_primes = 2*int(sys.argv[3])   #accounts for i%2 being skipped
    
    stop_primes = len(all_primes)
    if len(sys.argv) > 4:
        if int(sys.argv[4]) > 0:
            stop_primes = 2*int(sys.argv[4])    #accounts for i%2 being skipped

    logging = False
    if len(sys.argv) > 5:
        if sys.argv[5] == "log":
            logging = True

    assert stop_primes > start_primes

    primes = all_primes[start_primes:stop_primes]

    if logging:
        print(f'logging is turned ON, results are written to file results_ver{variant}_sam{num_samples}.txt\n\n')
    else:
        print("logging is turned OFF, results are not written to file\n\n")

    file = f'results_ver{variant}_sam{num_samples}.txt'

    ca = cost_graph(variant, num_samples, primes, file)
    ca2 = sum_f2(ca)






            