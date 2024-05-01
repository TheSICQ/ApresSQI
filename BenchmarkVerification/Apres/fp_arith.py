from random import randint

####################################################################################################################
################################## This file contains the functions needed to do  ##################################
##################################   the Fp and Fp^2 arithmetic (with counters)   ##################################
##################################  Also contains, and can change global params   ##################################
####################################################################################################################


#####################
# Fp arithmetic
#####################


def fp_add(a,b):
	"""
	Function to add two Fp elements a,b
	
	Input: integers a,b (mod p)
	Output: a+b (mod p)
	"""

	counter_fp[2] += 1

	return (a+b) % p

def fp_sub(a,b):
	"""
	Function to subtract two Fp elements a,b
	
	Input: integers a,b (mod p)
	Output: a-b (mod p)
	"""

	counter_fp[3] += 1

	return (a-b) % p

def fp_mul(a,b):
	"""
	Function to multiply two Fp elements a,b
	
	Input: integers a,b (mod p)
	Output: a*b (mod p)
	"""

	if (a == b):
		counter_fp[1] += 1
		return (a**2) % p

	counter_fp[0] += 1
	return (a*b) % p

def fp_sqr(a):
	"""
	Function to square an Fp element a
	
	Input: integers a (mod p)
	Output: a^2 (mod p)
	"""

	counter_fp[1] += 1

	return (a**2) % p

def fp_inv(a):
	"""
	Function to find the inverse of an Fp element a
	
	Input: integers a (mod p)
	Output: a^(-1) (mod p) such that a*a^(-1) = 1 (mod p)
	"""

	counter_fp[4] += 1

	return pow(a, p-2 , p)

def fp_neg(a):
	"""
	Function to find the negative of an Fp element a
	
	Input: integers a (mod p)
	Output: -a (mod p)
	"""

	return (-a) % p

def fp_issquare(a):
	"""
	Function to decide if an Fp element a is a square
	
	Input: integers a (mod p)
	Output: True if a is a square mod p, False otherwise
	"""

	counter_fp[5] += 1

	return pow(a, (p-1)//2, p) == 1

def fp_squareroot(a):
	"""
	Function to find the squareroot of an Fp element a
	
	Input: integers a (mod p)
	Output: square root b such that b^2 = a (mod p)
	"""

	counter_fp[6] += 1

	return pow(a, (p+1)//4, p)

def fp_expon(a,n):
	"""
	Function to compute the exponenetiation of an Fp element a by an integer
	
	Input: integer a (mod p) and integer n
	Output: a^n (mod p)
	"""
	
	nbits = [int(x) for x in bin(n)[2:]]    
	x = a

	for i in range(1, len(nbits)):  
		x = fp_sqr(x)
		if nbits[i] == 1:
			x = fp_mul(x,a)   
	return x


############################################
# Fp^2 arithmetic (i^2 + 1 = 0)
# We represent a1 + i*a2 \in Fp^2 as [a1,a2]
############################################


def fp2_add(a,b):
	"""
	Function to add two elements a,b in Fp2
	
	Input: pairs of integers a=[a1,a2],b=[b1,b2]
	Output: a+b in Fp^2
	"""

	counter_fp2[2] += 1

	return [fp_add(a[0],b[0]), fp_add(a[1],b[1])]

def fp2_sub(a,b):
	"""
	Function to subtract two elements a,b in Fp2
	
	Input: pairs of integers a=[a1,a2],b=[b1,b2]
	Output: a-b in Fp^2
	"""

	counter_fp2[3] += 1

	return [fp_sub(a[0],b[0]), fp_sub(a[1],b[1])]

def fp2_mul(a,b):
	"""
	Function to multiply two elements a,b in Fp2
	
	Input: pairs of integers a=[a1,a2],b=[b1,b2]
	Output: a*b in Fp^2
	"""

	counter_fp2[0] += 1

	t0 = fp_mul(a[0],b[0])
	t1 = fp_mul(a[1],b[1])
	t2 = fp_add(a[0],a[1])
	t3 = fp_add(b[0],b[1])

	res0 = fp_sub(t0,t1)
	t0 = fp_add(t0,t1)
	t1 = fp_mul(t2,t3)
	res1 = fp_sub(t1,t0)

	return [res0, res1]

def fp2_sqr(a):
	"""
	Function to square an element a in Fp2
	
	Input: pair of integers a=[a1,a2]
	Output: a^2 in Fp^2
	"""

	counter_fp2[1] += 1

	t0 = fp_add(a[0],a[1])
	t1 = fp_sub(a[0],a[1])
	t2 = fp_add(a[0],a[0])

	return [fp_mul(t0,t1), fp_mul(t2,a[1])]

def fp2_inv(a):
	"""
	Function to invert an element a in Fp2
	
	Input: pair of integers a=[a1,a2]
	Output: a^(-1) in Fp^2
	"""

	counter_fp2[4] += 1

	t0 = fp_sqr(a[0])
	t1 = fp_sqr(a[1])
	t0 = fp_add(t0,t1)
	t0 = fp_inv(t0)

	return([fp_mul(a[0],t0), fp_neg(fp_mul(a[1],t0))])

def fp2_neg(a):
	"""
	Function to negate an element a in Fp2
	
	Input: pair of integers a=[a1,a2]
	Output: -a in Fp^2
	"""

	return [fp_neg(a[0]), fp_neg(a[1])]

def fp2_issquare(a):
	"""
	Function to detect whether an element a is a square in Fp2
	
	Input: pair of integers a=[a1,a2]
	Output:  True if a is a square in Fp^2, False otherwise
	"""
	
	counter_fp2[5] += 1
	counter_fp[5] += 1

	t0 = fp_sqr(a[0])
	t1 = fp_sqr(a[1])
	t0 = fp_add(t0,t1)

	return pow(t0, (p-1)//2, p) == 1

def fp2_copy(a):
	"""
	Function that copies array that represents an Fp^2 element
	"""

	return [a[0], a[1]]


####################################################################
# Fast square roots in Fp2 
#
# Follows Scott's "A note on the calculation of some 
#       functions in finite fields: Tricks of the trade"
#
# Code mostly taken from https://github.com/microsoft/SuperSolver
#
####################################################################

"""
	MIT Licence

	Copyright (c) Microsoft Corporation. M. Corte-Real Santos, C. Costello, J. Shi

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE

"""

def tonelli_shanks_update(a,b,cs):
	"""
	Function that updates variables a,b according to Tonelli-Shanks Algorithm 
	
	Input: pairs of integers a=[a1,a2], b=[b1,b2] and Tonelli-Shanks constants cs
	Output: updated variables x,y
	"""

	
	c0 = int(cs[0])
	c1 = cs[1]
	c2 = cs[2]
	c3 = cs[3]
	x = a
	y = b

	for _ in range(2, c0 + 1):
		x = fp_sqr(x)
		y = fp_sqr(y)

	y = fp_sqr(y)
	y = fp_mul(y, x)

	if y != 1:
		a = fp_mul(a, c2)
		b = fp_mul(b, c3)

	x = fp_mul(a, b)
	z = fp_mul(x, b)

	for i in range(c0, 1, -1):
		w = z
		for _ in range(3, i+1):
			w = fp_sqr(w)
		if w != 1:
			x = fp_mul(x, c1)
		c1 = fp_sqr(c1)
		if w != 1:
			z = fp_mul(z, c1)
	return x, y


def fp2_squareroot(a):
	"""
	Function to find the squareroot of an element a in Fp2
	
	Input: pair of integers a=[a1,a2], representing a square in Fp^2 
	Output: square root b = [b1,b2] s.t. b^2 = a (mod p)
	"""
 
	counter_fp2[6] += 1

	# Get precomputed Tonelli-Shanks constants in globals
	c0 = constants[0]
	c1 = constants[1]
	c2 = constants[2]
	c3 = constants[3]
	c4 = constants[4]
	c5 = constants[5]
	c6 = constants[6]
	c7 = constants[7]
	c8 = constants[8]


	a1 = a[0]
	a2 = a[1]
	
	if a2 == 0:
		if fp_issquare(a1):
			return [fp_squareroot(a1), 0]
		else:
			a1_neg = fp_neg(a1)
			return [0, fp_squareroot(a1_neg)]

	t0 = a[0]
	t1 = a[1]

	t2 = fp_sqr(a1)
	t3 = fp_sqr(a2)
	t3 = fp_mul(t3, c0)
	t2 = fp_sub(t2, t3)

	t3 = fp_expon(t2, c4)

	t2, t3 = tonelli_shanks_update(t2, t3, [c2,c3,c0,c5]) 

	t3 = fp_add(a1, t2)
	t3 = fp_mul(t3, c1)

	t0 = fp_expon(t3, c4)

	t2 = t0

	for _ in range(1, c7 + 1):
		t2 = fp_sqr(t2)

	t2 = fp_mul(t2, c1)
	t2 = fp_mul(t2, t1)

	t1 = fp_expon(t3, c6)

	t2 = fp_mul(t1, t2)

	t3, t0 = tonelli_shanks_update(t3, t0, [c2,c3,c0,c5]) 

	t2 = fp_mul(t2, t3)

	if t0 == 1:
		t0 = t3
		t3 = t2
	else:
		t0 = t2
		t3 = fp_mul(t3, c8)

	b = [t0, t3]

	## Make the correct choice of squareroot (as per SQIsign specification)
	if b[0] % 2 != 0: 
		return [fp_mul(t0,c0), fp_mul(t3,c0)] 
	return b


################################################
# Get non-squares / squares functions
################################################

def random_nonsquare(i):
	return fp2_mul(fp2_mul(Z, [i,1]), [i,1])

def small_nonsquare(i):
	return small_nonsquares[i]

def small_square_with_minus_one_nonsquare(i):
	return small_square_min_one_nonsquares[i]

######################################################
# Other get-calls 
######################################################

def get_p():
	return p

def get_constants():
	return constants

##################
# GLOBALS
##################

from math import log

def log2(x):
	return log(x,2)

def update_p(newp):
	global counter_fp, counter_fp2, p, metric_fp, Z, small_nonsquares, small_square_min_one_nonsquares, constants, cost_isog, addition_chain, addition_chain_length, odd_factors
	cost_isog = 0
	counter_fp = [0, 0, 0, 0, 0, 0, 0]
	counter_fp2 = [0, 0, 0, 0, 0, 0, 0]
	p = newp
	Z = [randint(0, p), randint(0,p)]
	while True:
		Z = [randint(0, p), randint(0,p)]
		if fp2_issquare(Z):
			continue
		break

	small_nonsquares = []
	n = 0
	while len(small_nonsquares) < 2**8:
		if not fp2_issquare([n,1]):
			small_nonsquares.append([n,1])
		n += 1

	small_square_min_one_nonsquares = []
	n = 1
	while len(small_square_min_one_nonsquares) < 2**8:
		if fp2_issquare([n,1]) and not fp2_issquare([n-1,1]):
			small_square_min_one_nonsquares.append([n,1])
		n += 1
				#[mul, sqr, add, sub, inv, Legendre, sqrt]
	metric_fp = [1.0, 1.0, 0, 0, log2(p), log2(p), log2(p)]
	cost_isog = 0


	## Update the Tonelli-Shanks constants 

	c0 = -1 # tries c0 = -1,-2,... for quadratic extension
	while pow(c0, (p-1)//2, p) == 1:
		c0 -= 1
	c0 = c0 % p  # smallest QNR in GF(p)
	c1 = pow(2, p-2 , p) # constant for division by 2  
	c2 = 1 # Tonelli-Shanks "e" in Michael Scott's paper
	while (p-1) % pow(2, c2 + 1, p) == 0:
		c2 += 1
	c3 = pow(c0, (p-1)*pow(pow(2, c2, p), p-2, p), p) # "z" in Michael Scott's paper
	c4 = (p-1-pow(2, c2, p))*pow(pow(2, c2 + 1, p), p-2, p) # y = x^c4 in Michael Scott's paper
	c4 = int(c4 % p)
	c6 = pow(2, c2, p) - 1
	c7 = c2 + 1
	c8 = pow(c0,p-2,p) #1/c0 in Michael Scott's paper
	c5 = pow(c0, c4, p)
	constants = [c0,c1,c2,c3,c4,c5,c6,c7,c8]

	# pre-computed differential addition chains
	addition_chain = {149: 416, 257: 641, 131: 64, 223: 708, 5: 0, 263: 592, 103: 48, 137: 424, 139: 388, 13: 0, 271: 648, 17: 10, 19: 4, 277: 672, 151: 180, 281: 328, 283: 576, 157: 384, 31: 8, 163: 296, 503: 320, 293: 514, 167: 266, 41: 42, 43: 34, 173: 192, 307: 264, 47: 32, 179: 272, 53: 106, 311: 160, 23: 24, 313: 80, 347: 1384, 59: 88, 317: 320, 181: 258, 193: 130, 67: 72, 197: 785, 71: 80, 73: 20, 587: 2832, 239: 682, 269: 552, 11: 4, 337: 32, 211: 788, 7: 2, 89: 0, 331: 1704, 79: 16, 349: 1448, 199: 256, 37: 48, 3: 0, 353: 1424, 227: 776, 101: 164, 367: 1410, 359: 1600, 233: 0, 97: 192, 107: 96, 29: 16, 109: 132, 61: 81, 229: 784, 113: 464, 83: 210, 373: 1354, 251: 554, 241: 832, 191: 136, 127: 417}
	addition_chain_length = {257: 10, 131: 8, 17: 4, 179: 9, 263: 10, 137: 9, 139: 9, 13: 3, 269: 10, 271: 10, 43: 6, 107: 8, 277: 10, 151: 9, 587: 12, 281: 10, 283: 10, 5: 1, 157: 9, 31: 5, 163: 9, 293: 10, 167: 9, 41: 6, 71: 7, 353: 11, 173: 9, 47: 6, 359: 11, 53: 7, 307: 10, 181: 9, 311: 10, 23: 5, 313: 10, 59: 7, 317: 10, 191: 9, 19: 4, 193: 9, 67: 7, 197: 10, 199: 9, 73: 7, 331: 11, 79: 7, 11: 3, 337: 10, 37: 6, 211: 10, 7: 2, 89: 7, 347: 11, 101: 8, 349: 11, 223: 10, 3: 0, 97: 8, 61: 7, 227: 10, 113: 9, 229: 10, 367: 11, 103: 8, 233: 9, 127: 9, 29: 5, 109: 8, 239: 10, 241: 10, 83: 8, 373: 11, 251: 10, 149: 9, 503: 11}
 
	reset_counter()
	return

def reset_counter():
	global counter_fp, counter_fp2, cost_isog
	cost_isog = 0
	counter_fp = [0, 0, 0, 0, 0, 0, 0]
	counter_fp2 = [0, 0, 0, 0, 0, 0, 0]
	return

def get_cost():
	return sum([op*cost for op, cost in zip(counter_fp, metric_fp)])

def get_addition_chain_and_length(l):
    return addition_chain[l], addition_chain_length[l]

def get_addition_chain_primes():
    return addition_chain.keys()