
####################################################################################################################
###############################   This file contains an example to run  (see README)  ##############################
####################################################################################################################


from ApresSQI import SQIsign
import time
from math import ceil

from sys import argv

def compute_size(sigma):
    def flatten_list(input_list):
        flattened = []
        for item in input_list:
            if isinstance(item, list) or isinstance(item, tuple):
                flattened.extend(flatten_list(item))
            else:
                flattened.append(item)
        return flattened
    if len(sigma) == 3: #unseeded
        seeds = []
    elif len(sigma) == 4: # uncompressed
        seeds = sigma[-1] # Hash seeds
        sigma = sigma[:-1]
    elif len(sigma) == 5: #seeded
        seeds = sigma[-2:] #both seeds for zip and chall
        sigma = sigma[:-2]
    nums = [int(n) for n in flatten_list(sigma)]
    bit_len = sum([(1 if n == 0 else n.bit_length()) for n in nums])
    bit_len += len(flatten_list(seeds))*8 # One byte per seed
    return ceil(bit_len/8)


if __name__=="__main__":
    param = 'toy'
    if len(argv) > 1:
        param = argv[1]
    seeded = True
    if len(argv) > 2:
        seeded = not argv[2] == 'False'
    compressed = True
    if len(argv) > 3:
        compressed = not argv[3] == 'False'

    assert param in ['toy', 'NIST', '7-block', '4-block']

    msg = "Ska'ru værra med i bakken?"
    
    topstr = f'############# Using {param} prime #############'
    print(topstr)

    signer = SQIsign(param)
    signer.load_privkey(param=param)

    print(f'Signing message: \n{msg}')
    print('\nSignature will be ' + ('' if compressed else 'un') +'compressed', end='')
    if compressed:
        print(' and ' + ('' if seeded else 'un') +'seeded')
    else:
        print()
    print('#'*len(topstr) + '\n\n')

    tstart = time.time()
    sigma = signer.Sign(msg, seeded = seeded, compressed = compressed)

    total_time = time.time() - tstart


    print('\n\n' + '#'*len(topstr) + '\n')
    print(f'Done! signature was {sigma}' + '\n')

    print("Verifying...")
    verified = signer.verify(msg, sigma, signer.pk)
    print('Signature was ' + ('CORRECT' if verified else 'WRONG'))

    size_str = f"""
    ################################################
    #
    #           signature size: {compute_size(sigma)} B
    #
    ################################################n
    """

    time_str = f"""
    ################################################
    #
    #         signing time: {round(total_time,3)} seconds
    #
    ################################################
    """
    print(size_str)
    print(time_str)






