from random import randint

####################################################################################################################
################################## This file contains the functions needed to do  ##################################
##################################   the Fp and Fp^2 arithmetic (with counters)   ##################################
####################################################################################################################

####################################################################################################################
###########################$$###### Also contains, and can change global params  ###################################
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

def fp2_squareroot(a):
	"""
	Function to find the squareroot of an element a in Fp2
	
	Input: pair of integers a=[a1,a2], representing a square in Fp^2 
	Output: square root b = [b1,b2] s.t. b^2 = a (mod p)
	"""
	
	counter_fp2[6] += 1

	if a[1] == 0:
		if fp_issquare(a[0]):
			return [fp_squareroot(a[0]), 0]
		else:
			t0 = fp_neg(a[0])
			return [0, fp_squareroot(t0)]

	sdelta = fp_sqr(a[0])
	t1 = fp_sqr(a[1])
	sdelta = fp_add(sdelta, t1)
	sdelta = fp_squareroot(sdelta)

	inv2 = fp_inv(2)
	re = fp_add(a[0], sdelta)
	re = fp_mul(re, inv2)

	if not fp_issquare(re):
		re = fp_sub(a[0], sdelta)
		re = fp_mul(re, inv2)

	re = fp_squareroot(re)
	im = fp_inv(re)
	im = fp_mul(im, inv2)
	im = fp_mul(im, a[1])

	return [re, im]


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

##################
# GLOBALS
##################

from math import log

def log2(x):
	return log(x,2)

def update_p(newp):
	global counter_fp, counter_fp2, p, metric_fp, Z, small_nonsquares, small_square_min_one_nonsquares, cost_isog
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

def reset_counter():
	global counter_fp, counter_fp2, cost_isog
	cost_isog = 0
	counter_fp = [0, 0, 0, 0, 0, 0, 0]
	counter_fp2 = [0, 0, 0, 0, 0, 0, 0]
	return

def get_cost():
	return sum([op*cost for op, cost in zip(counter_fp, metric_fp)])