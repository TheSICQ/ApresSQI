from ApresSQI import SQIsign

import time

def time_sign(param):
    signer = SQIsign(param)
    signer.load_privkey(param=param)
    msg = "Øl, øl, øl og aftersqi!"
    t_start = time.time()
    sigma = signer.Sign(msg)
    tot_time = time.time() - t_start
    print(f"Signing finished: Took a total of {tot_time} seconds")
    print(sigma)
    ("Verifying...")
    verified = signer.verify(msg, sigma, signer.pk)
    print('Signature was ' + ('CORRECT' if verified else 'WRONG'))
    assert verified
    return tot_time

if __name__ == "__main__":
    #params = ['NIST', '7-block', '4-block']
    params = ['NIST', '8513034219037441780170691209753296498696014329521974009944792576819199999999']
    timings = {}
    reps = 3
    for p in params:
        timings[p] = []
        for _ in range(reps):
            timings[p].append(time_sign(p))
        print("Data now:")
        print(timings)
    for p, t in timings.items():
        print(f"param: {p}:")
        print(f"average time: {round(sum(t)/len(t), 5)}\n\n")