from isogenychains import fp_issquare, fp_inv 

def generate_constants(p):
    # inputs: (Global) p
    # outputs: Tonelli-Shanks constants for p
    c0 = -1 # tries c0 = -1,-2,... for quadratic extension
    while fp_issquare(c0):
        c0 -= 1
	

    c0 = c0 % p  # smallest QNR in GF(p)
    c1 = fp_inv(2) # constant for division by 2  
    c2 = 1 # Tonelli-Shanks "e" in Michael Scott's paper
    while (p-1) % pow(2, c2 + 1, p) == 0:
        c2 += 1
	
     
    c3 = pow(c0, (p-1)*fp_inv(pow(2, c2, p)), p) # "z" in Michael Scott's paper
    c4 = (p-1-pow(2, c2, p))*fp_inv(pow(2, c2 + 1, p)) # y = x^c4 in Michael Scott's paper
    c4 = int(c4 % p)
    c6 = pow(2, c2, p) - 1
    c7 = c2 + 1
    c8 = fp_inv(c0) #1/c0 in Michael Scott's paper
    c5 = pow(c0, c4, p)
    return [c0,c1,c2,c3,c4,c5,c6,c7,c8]