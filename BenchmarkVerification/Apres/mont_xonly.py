from fp_arith import fp2_sub, fp2_add, fp2_copy, fp2_inv, fp2_issquare, fp2_mul, fp2_neg, fp2_sqr, fp2_squareroot, fp_add, fp_expon, fp_inv, fp_issquare, fp_mul, fp_neg, fp_sqr, fp_squareroot, fp_sub, update_p, reset_counter, small_nonsquare, small_square_with_minus_one_nonsquare, random_nonsquare, get_addition_chain_and_length, get_addition_chain_primes

from testing import valuation
####################################################################################################################
################################## This file contains the functions needed to do  ##################################
##################################    x-only Montgomery Scalar Multiplication     ##################################
####################################################################################################################

##################
# DOUBLING
##################


def xDBL(P, ProjA):
	"""
	This function does Montgomery doubling with projective curve constant (A:C)
	
	Input: projective coordinates xP=P[0]/P[1]
		   curve constant A = ProjA[0]/ProjA[1]
	Output: projective coordinates x(2P)=X2/Z2
	"""

	t0 = fp2_sub(P[0], P[1])
	t1 = fp2_add(P[0], P[1])
	t0 = fp2_sqr(t0)
	t1 = fp2_sqr(t1)
	Z2 = fp2_mul(ProjA[1], t0)
	Z2 = fp2_add(Z2, Z2)
	Z2 = fp2_add(Z2, Z2)
	X2 = fp2_mul(Z2, t1)
	t1 = fp2_sub(t1, t0)
	t0 = fp2_add(ProjA[1], ProjA[1])
	t0 = fp2_add(t0, ProjA[0])
	t0 = fp2_mul(t0, t1)
	Z2 = fp2_add(Z2, t0)
	Z2 = fp2_mul(Z2, t1)

	return [X2, Z2]   #cost: 4M+2S+8a

def xDBLaff(P, A):
	"""
	This function does Montgomery doubling with affine curve constant
	
	Input: projective coordinates xP=P[0]/P[1]
		   curve constant A = ProjA[0]/1
	Output: projective coordinates x(2P)=X2/Z2
	"""

	assert A[1] == [1,0]

	t0 = fp2_sub(P[0], P[1])
	t1 = fp2_add(P[0], P[1])
	t0 = fp2_sqr(t0)
	t1 = fp2_sqr(t1)
	#Z2 = fp2_mul(ProjA[1], t0)
	Z2 = fp2_add(t0, t0)
	Z2 = fp2_add(Z2, Z2)
	X2 = fp2_mul(Z2, t1)
	t1 = fp2_sub(t1, t0)
	t0 = fp2_add([1,0], [1,0])
	t0 = fp2_add(t0, A[0])
	t0 = fp2_mul(t0, t1)
	Z2 = fp2_add(Z2, t0)
	Z2 = fp2_mul(Z2, t1)

	return [X2, Z2]   #cost: 3M+2S+8a


def xDBLA24(P, A24):
	"""
	Compute Montgomery doubling with projective curve constant (A+2C : 4C)
	Follows code from SQIsign NIST version 1.0 (xDBLv2)
	
	Input: projective coordinates xP=P[0]/P[1]
		   curve constant A24plus = (A+2C : 4C)
	Output: projective coordinates x(2P)=X2/Z2
	"""

	t0 = fp2_add(P[0], P[1])
	t0 = fp2_sqr(t0)
	t1 = fp2_sub(P[0], P[1])
	t1 = fp2_sqr(t1)
	t2 = fp2_sub(t0, t1)
	t1 = fp2_mul(t1, A24[1])
	X2 = fp2_mul(t0, t1)
	t0 = fp2_mul(t2, A24[0])
	t0 = fp2_add(t0, t1)
	Z2 = fp2_mul(t0, t2)

	return [X2, Z2]   #cost: 4M+2S+4a


def xDBLe(P, ProjA, e):
	"""
	Comptutes e Montgomery doublings using xDBL
	
	Input: projective coordinates xP=P[0]/P[1]
		   curve constant A = ProjA[0]/ProjA[1]
	Output: projective coordinates x([2^e]P)=X2/Z2
	"""

	Pout = fp2_copy(P)

	for _ in range(e):
		Pout = xDBL(Pout,ProjA)

	return Pout

def xDBLe_aff(P, A, e):
	"""
	Comptutes e Montgomery doublings using xDBLaff
	
	Input: projective coordinates xP=P[0]/P[1]
		   curve constant A = ProjA[0]/1
	Output: projective coordinates x([2^e]P)=X2/Z2
	"""

	assert A[1] == [1,0]

	Pout = fp2_copy(P)

	for _ in range(e):
		Pout = xDBLaff(Pout,A)

	return Pout

def xDBLeA24(P, ProjA24, e):
	"""
	Comptutes e Montgomery doublings using xDBLA24
	
	Input: projectivee coordinates xP=P[0]/P[1]
		   curve constant A = ProjA[0]/ProjA[1]
	Output: projective coordinates x([2^e]P)=X2/Z2
	"""

	Pout = fp2_copy(P)

	for _ in range(e):
		Pout = xDBLA24(Pout, ProjA24)

	return Pout

##################
# TRIPLING
##################

def xTPL(P, ProjA):
	"""
	Computes Montgomery tripling with projective curve constant
	Code from the latest SIKE spec
	
	Input: projective coordinates xP=P[0]/P[1]
		   curve constant A = ProjA[0]/ProjA[1]
	Output: projective coordinates x(2P)=X2/Z2
	"""
	
	C3 = fp2_add(ProjA[1], ProjA[1])
	A3 = fp2_add(ProjA[0], C3)
	C3 = fp2_sub(ProjA[0], C3)

	t0 = fp2_sub(P[0], P[1])
	t2 = fp2_sqr(t0)
	t1 = fp2_add(P[0], P[1])
	t3 = fp2_sqr(t1)
	t4 = fp2_add(t1, t0)
	t0 = fp2_sub(t1, t0)
	t1 = fp2_sqr(t4)
	t1 = fp2_sub(t1, t3)
	t1 = fp2_sub(t1, t2)
	t5 = fp2_mul(t3, A3)
	t3 = fp2_mul(t3, t5)
	t6 = fp2_mul(t2, C3)
	t2 = fp2_mul(t2, t6)
	t3 = fp2_sub(t2, t3)
	t2 = fp2_sub(t5, t6)
	t1 = fp2_mul(t2, t1)
	t2 = fp2_add(t3, t1)
	t2 = fp2_sqr(t2)
	X3 = fp2_mul(t2, t4)
	t1 = fp2_sub(t3, t1)
	t1 = fp2_sqr(t1)
	Z3 = fp2_mul(t1, t0)

	return [X3, Z3]

def xTPLA24pm(P, A24pm):
	"""
	Computes Montgomery tripling with projective curve constant A24plusminus
	Code from latest SIKE spec
	
	Input: projective coordinates xP=P[0]/P[1]
		   curve constant A24pm = (A+2C : A-2C)
	Output: projective coordinates x(2P)=X2/Z2
	"""
	
	t0 = fp2_sub(P[0], P[1])
	t2 = fp2_sqr(t0)
	t1 = fp2_add(P[0], P[1])
	t3 = fp2_sqr(t1)
	t4 = fp2_add(t1, t0)
	t0 = fp2_sub(t1, t0)
	t1 = fp2_sqr(t4)
	t1 = fp2_sub(t1, t3)
	t1 = fp2_sub(t1, t2)
	t5 = fp2_mul(t3, A24pm[0])
	t3 = fp2_mul(t3, t5)
	t6 = fp2_mul(t2, A24pm[1])
	t2 = fp2_mul(t2, t6)
	t3 = fp2_sub(t2, t3)
	t2 = fp2_sub(t5, t6)
	t1 = fp2_mul(t2, t1)
	t2 = fp2_add(t3, t1)
	t2 = fp2_sqr(t2)
	X3 = fp2_mul(t2, t4)
	t1 = fp2_sub(t3, t1)
	t1 = fp2_sqr(t1)
	Z3 = fp2_mul(t1, t0)

	return [X3, Z3]

def xTPLe(P, ProjA, e):
	"""
	Computes e Montgomery triplings using xTPL
	Code from latest SIKE spec
	
	Input: projective coordinates xP=P[0]/P[1]
		   curve constant A = ProjA[0]/ProjA[1]
	Output: projective coordinates x(2P)=X2/Z2
	"""
	Pout = fp2_copy(P)

	for _ in range(e):
		Pout = xTPL(Pout,ProjA)

	return Pout

def xTPLeA24pm(P, A24pm, e):
	"""
	Computes e Montgomery triplings using xTPLA24pm
	Code from latest SIKE spec
	
	Input: projective coordinates xP=P[0]/P[1]
		   curve constant A24pm = (A+2C : A-2C)
	Output: projective coordinates x(2P)=X2/Z2
	"""

	Pout = fp2_copy(P)

	for _ in range(e):
		Pout = xTPLA24pm(Pout, A24pm)

	return Pout


########################
# DIFFERENTIAL ADDITION
########################

def xADD(P, Q, PmQ):
	"""
	Computes Montgomery differential addition
	
	Input: projective points P, Q
		   projective difference PmQ = P-Q
	Output: projective coordinates of P+Q
	"""
 
	t0 = fp2_add(P[0], P[1])
	t1 = fp2_sub(P[0], P[1])
	t2 = fp2_sub(Q[0], Q[1])
	t3 = fp2_add(Q[0], Q[1])
	t0 = fp2_mul(t2, t0)
	t1 = fp2_mul(t3, t1)
	t3 = fp2_sub(t0, t1)
	t2 = fp2_add(t0, t1)
	t3 = fp2_sqr(t3)
	XQP = fp2_sqr(t2)
	ZQP = fp2_mul(PmQ[0], t3)
	XQP = fp2_mul(XQP, PmQ[1])

	return [XQP, ZQP]    #cost: 4M+2S+6a

def small_mul(a, b):

	"""
		Multiply a,b where a = [x,0] with small x
		Input: a = [x,0] with small x, b in Fp^2
		Output: a*b
	"""

	assert a[0] < 128
	assert a[1] <= 1

	res = b.copy()
	for _ in range(a[0]-1):
		res = fp2_add(res, b)	

	if a[1] == 0:
		return res
	elif a[1] == 1:
		res = fp2_add(res, [fp_neg(b[1]), b[0]])
		return res
	else:
		assert False


def xADDsmall(P, Q, PmQ):
	"""
	Input: projective points P, Q (with small TODO?)
		   projective difference PmQ = P-Q with PmQ[1] = [1,0] (i.e., affine)
	Output: projective coordinates of P+Q
	"""

	assert PmQ[1] == [1,0]

	t0 = fp2_add(P[0], P[1])
	t1 = fp2_sub(P[0], P[1])
	t2 = fp2_sub(Q[0], Q[1])
	t3 = fp2_add(Q[0], Q[1])
	t0 = fp2_mul(t2, t0)
	t1 = fp2_mul(t3, t1)
	t3 = fp2_sub(t0, t1)
	t2 = fp2_add(t0, t1)
	t3 = fp2_sqr(t3)
	XQP = fp2_sqr(t2)

	ZQP = small_mul(PmQ[0], t3)

	#assert ZQP == fp2_mul(PmQ[0], t3)		
	#XQP = fp2_mul(XQP, PmQ[1])

	return [XQP, ZQP]    #cost: 2M+2S+6A + cost of small_mul


def xDBLADD(P, Q, PmQ, ProjA):
	"""
	Function for step in Montgomery ladder
	Simultaneous doubling and differential addition
	
	Input: projective points P, Q
		   projective difference PmQ = P-Q
		   curve constant A = ProjA[0]/ProjA[1]
	Output: projective coordinates of [2]P and P+Q
	"""
	
	t0 = fp2_add(P[0], P[1])
	t1 = fp2_sub(P[0], P[1])
	X2P = fp2_sqr(t0)
	t2 = fp2_sub(Q[0], Q[1])
	XQP = fp2_add(Q[0], Q[1])
	t0 = fp2_mul(t0, t2)
	Z2P = fp2_sqr(t1)
	t1 = fp2_mul(t1, XQP)
	t2 = fp2_sub(X2P, Z2P)
	t3 = fp2_add(ProjA[1], ProjA[1])
	Z2P = fp2_mul(t3, Z2P)
	t3 = fp2_add(ProjA[0], t3)
	Z2P = fp2_add(Z2P, Z2P)
	X2P = fp2_mul(X2P, Z2P)
	XQP = fp2_mul(t3, t2)
	ZQP = fp2_sub(t0, t1)
	Z2P = fp2_add(XQP, Z2P)
	XQP = fp2_add(t0, t1)
	Z2P = fp2_mul(Z2P, t2)
	ZQP = fp2_sqr(ZQP)
	XQP = fp2_sqr(XQP)
	ZQP = fp2_mul(PmQ[0], ZQP)
	XQP = fp2_mul(XQP, PmQ[1])

	return [X2P, Z2P], [XQP, ZQP] # 4S+8M+11A

def xDBLADDaff(P, Q, PmQ, A):
	"""
	Function for step in Montgomery ladder with affine curve constant and P-Q
	Simultaneous doubling and differential addition
	
	Input: projective points P, Q
		   difference PmQ = P-Q with PmQ[1] = [1,0] (i.e., affine)
		   curve constant A = ProjA[0]/1
	Output: projective coordinates of [2]P and P+Q
	"""

	assert A[1] == [1,0]
	assert PmQ[1] == [1,0]

	t0 = fp2_add(P[0], P[1])
	t1 = fp2_sub(P[0], P[1])
	X2P = fp2_sqr(t0)
	t2 = fp2_sub(Q[0], Q[1])
	XQP = fp2_add(Q[0], Q[1])
	t0 = fp2_mul(t0, t2)
	Z2P = fp2_sqr(t1)
	t1 = fp2_mul(t1, XQP)
	t2 = fp2_sub(X2P, Z2P)
	#t3 = fp2_add(ProjA[1], ProjA[1]) this is just 2
	t3 = fp2_add([1,0], [1,0])
	#Z2P = fp2_mul(t3, Z2P)
	Z2P = fp2_add(Z2P, Z2P)
	t3 = fp2_add(A[0], t3)
	Z2P = fp2_add(Z2P, Z2P)
	X2P = fp2_mul(X2P, Z2P)
	XQP = fp2_mul(t3, t2)
	ZQP = fp2_sub(t0, t1)
	Z2P = fp2_add(XQP, Z2P)
	XQP = fp2_add(t0, t1)
	Z2P = fp2_mul(Z2P, t2)
	ZQP = fp2_sqr(ZQP)
	XQP = fp2_sqr(XQP)
	ZQP = fp2_mul(PmQ[0], ZQP)
	# XQP = fp2_mul(XQP, PmQ[1])

	return [X2P, Z2P], [XQP, ZQP] # 4S+7M+11A


def xDBLADDaffs(P, Q, PmQ, A):
	"""
	Function for step in Montgomery ladder with affine curve constant
	Simultaneous doubling and differential addition
	
	Input: projective points P, Q
		   projective difference PmQ = P-Q
		   curve constant A = ProjA[0]/1
	Output: projective coordinates of [2]P and P+Q
	"""

	assert A[1] == [1,0]

	t0 = fp2_add(P[0], P[1])
	t1 = fp2_sub(P[0], P[1])
	X2P = fp2_sqr(t0)
	t2 = fp2_sub(Q[0], Q[1])
	XQP = fp2_add(Q[0], Q[1])
	t0 = fp2_mul(t0, t2)
	Z2P = fp2_sqr(t1)
	t1 = fp2_mul(t1, XQP)
	t2 = fp2_sub(X2P, Z2P)
	#t3 = fp2_add(ProjA[1], ProjA[1]) this is just 2
	t3 = fp2_add([1,0], [1,0])
	#Z2P = fp2_mul(t3, Z2P)
	Z2P = fp2_add(Z2P, Z2P)
	t3 = fp2_add(A[0], t3)
	Z2P = fp2_add(Z2P, Z2P)
	X2P = fp2_mul(X2P, Z2P)
	XQP = fp2_mul(t3, t2)
	ZQP = fp2_sub(t0, t1)
	Z2P = fp2_add(XQP, Z2P)
	XQP = fp2_add(t0, t1)
	Z2P = fp2_mul(Z2P, t2)
	ZQP = fp2_sqr(ZQP)
	XQP = fp2_sqr(XQP)
	ZQP = fp2_mul(PmQ[0], ZQP)
	XQP = fp2_mul(XQP, PmQ[1])

	return [X2P, Z2P], [XQP, ZQP] # 4S+7M+11A

def xDBLADDsmallP(P, Q, PmQ, A):
	"""
	Function for step in Montgomery ladder with affine curve constant and Q and affine, small P
	Simultaneous doubling and differential addition
	
	projective point P, with P[1] = [1,0] (i.e., affine)
			scalar m, curve constant A = ProjA[0]/1
			
	Input: projective points P, Q (with small P, and P[1] = Q[1] = [1,0], i.e., affine)
		   projective difference PmQ = P-Q 
		   curve constant A = ProjA[0]/1
	Output: projective coordinates of [2]P and P+Q
	"""

	assert P[0][0] <= 128
	assert P[0][1] <= 128
	assert P[1] == [1,0]
	assert Q[1] == [1,0]
	assert A[1] == [1,0]

	t0 = fp2_add(P[0], P[1])
	assert t0 == P[0]
	t1 = fp2_sub(P[0], P[1])
	assert t0 == t1

	X2P = fp2_sqr(t0)
	t2 = fp2_sub(Q[0], Q[1])
	XQP = fp2_add(Q[0], Q[1])
	t0 = fp2_mul(t0, t2)
	Z2P = fp2_sqr(t1)
	t1 = fp2_mul(t1, XQP)
	t2 = fp2_sub(X2P, Z2P)
	#t3 = fp2_add(ProjA[1], ProjA[1]) this is just 2
	t3 = fp2_add([1,0], [1,0])
	#Z2P = fp2_mul(t3, Z2P)
	Z2P = fp2_add(Z2P, Z2P)
	t3 = fp2_add(A[0], t3)
	Z2P = fp2_add(Z2P, Z2P)
	X2P = fp2_mul(X2P, Z2P)
	XQP = fp2_mul(t3, t2)
	ZQP = fp2_sub(t0, t1)
	Z2P = fp2_add(XQP, Z2P)
	XQP = fp2_add(t0, t1)
	Z2P = fp2_mul(Z2P, t2)
	ZQP = fp2_sqr(ZQP)
	XQP = fp2_sqr(XQP)
	ZQP = fp2_mul(PmQ[0], ZQP)
	XQP = fp2_mul(XQP, PmQ[1])

	return [X2P, Z2P], [XQP, ZQP] # 4S+7M+11A

#####################
# MONTGOMERY LADDER
#####################

def xMUL_noDACs(P, m, ProjA):
	"""
	Computes x-only Montgomery Scalar Multiplication
	
	Input: projective point P,
			scalar m, curve constant A = ProjA[0]/ProjA[1]
	Output: projective coordinates of [m]P
	"""

	binary = lambda n: n>0 and [n&1]+binary(n>>1) or []
	bits = binary(m)

	P0 = [[1, 0], [0, 0]]
	P1 = [P[0], P[1]]

	for i in range(len(bits)-1, -1, -1):
		if bits[i] == 0:
			P0, P1 = xDBLADD(P0, P1, P, ProjA)
		else:
			P1, P0 = xDBLADD(P1, P0, P, ProjA)

	return P0


def xMULaff(P, m, A):
	"""
	Computes x-only Montgomery Scalar Multiplication with affine curve constant and P
	Input: projective point P, with P[1] = [1,0] (i.e., affine)
			scalar m, curve constant A = ProjA[0]/1
	Output: projective coordinates of [m]P
	"""

	binary = lambda n: n>0 and [n&1]+binary(n>>1) or []
	bits = binary(m)

	assert A[1] == [1,0]
	assert P[1] == [1,0]

	P0 = [[1, 0], [0, 0]]
	P1 = [P[0], P[1]]

	for i in range(len(bits)-1, -1, -1):
		if bits[i] == 0:
			P0, P1 = xDBLADDaff(P0, P1, P, A)
		else:
			P1, P0 = xDBLADDaff(P1, P0, P, A)

	return P0


def xMULaffs_noDACs(P, m, A):
	"""
	Computes x-only Montgomery Scalar Multiplication with affine curve constant
	Input: projective point P,
			scalar m, curve constant A = ProjA[0]/1
	Output: projective coordinates of [m]P
	"""


	binary = lambda n: n>0 and [n&1]+binary(n>>1) or []
	bits = binary(m)

	assert A[1] == [1,0]

	P0 = [[1, 0], [0, 0]]
	P1 = [P[0], P[1]]

	for i in range(len(bits)-1, -1, -1):
		if bits[i] == 0:
			P0, P1 = xDBLADDaffs(P0, P1, P, A)
		else:
			P1, P0 = xDBLADDaffs(P1, P0, P, A)

	return P0


def xMUL_prime(P, l, ProjA):
	"""
	Computes x-only Montgomery Scalar Multiplication for a prime scalar for which we have pre-computed the shortest differential addition chain
	
	Input: projective point P,
			prime scalar l, curve constant A = ProjA[0]/ProjA[1]
	Output: projective coordinates of [l]P
	"""

	# Get addition chain and its length from globals
	chain, chain_len = get_addition_chain_and_length(l)
 
 	#  Initial 3-tuple of points
	R0 = P
	R1 = xDBL(P, ProjA)
	R2 = xADD(R1, R0, P)
	R = [R0,R1,R2]

 	# Main loop
	for _ in range(chain_len):
		# if isinfinity(R[(chain & 0x1)]) == 1:
		# 	T = xDBL(R[2], ProjA)
		# else:
		T = xADD(R[2], R[(chain & 1) ^ 1], R[chain & 1])
		R = [R[(chain & 1) ^ 1], R[2], T]
		chain >>= 1

	lP = R[2]
	return lP


def xMUL(P, m, ProjA):
	"""
 	Computes x-only Montgomery Scalar Multiplication using pre-computed differential addition chains (DACs)
	
	Input: projective point P, scalar m, m_facs the odd factors of m, 
			f the maximum power of 2 dividing m,
			curve constant A = ProjA[0]/ProjA[1]
	Output: projective coordinates of [m]P
	"""

 	# for the primes that we have precomputed DACs, see which divide m
	ls = get_addition_chain_primes()
	facs = [l for l in ls if valuation(m,l) > 0]
	
	mP = P

	for l in facs:
		m = m//l
		mP = xMUL_prime(mP, l, ProjA)

	# get the power of 2 dividing m
	e = valuation(m, 2)
	if e > 0:
		m = m//2**e
		mP = xDBLe(mP, ProjA, e)

	# for the remaining factors do the normal xMUL
	mP = xMUL_noDACs(mP, m, ProjA)
	return mP




def xMULaffs_prime(P, l, ProjA):
	"""
	Computes x-only Montgomery Scalar Multiplication for a prime scalar for which we have pre-computed the shortest differential addition chain
	
	Input: projective point P,
			prime scalar l, curve constant A = ProjA[0]/1 
	Output: projective coordinates of [l]P
	"""

	# Get addition chain and its length from globals
	chain, chain_len = get_addition_chain_and_length(l)
 
 	#  Initial 3-tuple of points
	R0 = P
	R1 = xDBLaff(P, ProjA) #A is affine, P is projective so we use xDBLaff
	R2 = xADD(R1, R0, P) # affine A does not change xADD
	R = [R0,R1,R2]

 	# Main loop
	for _ in range(chain_len):
		# if isinfinity(R[(chain & 0x1)]) == 1:
		# 	T = xDBL(R[2], ProjA)
		# else:
		T = xADD(R[2], R[(chain & 1) ^ 1], R[chain & 1])
		R = [R[(chain & 1) ^ 1], R[2], T]
		chain >>= 1

	lP = R[2]
	return lP


def xMULaffs(P, m, ProjA):
	"""
 	Computes x-only Montgomery Scalar Multiplication using pre-computed differential addition chains (DACs)
	
	Input: projective point P, scalar m, m_facs the odd factors of m, 
			f the maximum power of 2 dividing m,
			curve constant A = ProjA[0]/1
	Output: projective coordinates of [m]P
	"""

 	# for the primes that we have precomputed DACs, see which divide m
	ls = get_addition_chain_primes()
	facs = [l for l in ls if valuation(m,l) > 0]
	
	mP = P

	for l in facs:
		m = m//l
		mP = xMULaffs_prime(mP, l, ProjA)

	# get the power of 2 dividing m
	e = valuation(m, 2)
	if e > 0:
		m = m//2**e
		mP = xDBLe_aff(mP, ProjA, e)

	# for the remaining factors do the normal xMUL
	mP = xMULaffs_noDACs(mP, m, ProjA)
	return mP



def xMULc(P, m, A):
	"""
	Computes x-only Montgomery Scalar Multiplication with affine curve constant and affine, small P
	Input: affine point P with small coordinates, with P[1] = [1,0]
			scalar m, curve constant A = ProjA[0]/1
	Output: projective coordinates of [m]P
	"""
	
	binary = lambda n: n>0 and [n&1]+binary(n>>1) or []
	bits = binary(m)

	assert A[1] == [1,0]
	assert P[1] == [1,0]
	assert P[0][0] < 128
	assert P[0][1] <= 1

	P0 = [[1, 0], [0, 0]]
	P1 = [P[0], P[1]]

	for i in range(len(bits)-1, -1, -1):
		if bits[i] == 0:
			P1 = xADDsmall(P0, P1, P)
			P0 = xDBLaff(P0, A)

		else:
			P0 = xADDsmall(P0, P1, P)
			P1 = xDBLaff(P1, A)

	return P0

#################
# BIDIM MULT
#################

def xMULbidim(P, Q, PmQ, m, n, ProjA):
	"""
	Computes Bidimensional Multiplication
	Algorithm 9 from eprint 2017/212

	Input: projective points P, Q, PmQ = P-Q
			scalars m, n
	Output: projective point R = [m]P + [n]Q
	"""

	# Following notation from Alg. 9
	s0, s1 = m, n
	x0, x1, xmin = fp2_copy(P), fp2_copy(Q), fp2_copy(PmQ)

	while s0 != 0:

		if s1 < s0:
			s0, s1 = s1, s0
			x0, x1 = x1, x0

		if s1 <= 4*s0:
			s1 = s1 - s0
			tmp = xADD(x1, x0, xmin)
			xmin = fp2_copy(x0)
			x0 = fp2_copy(tmp)
		elif s0 % 2 == s1 % 2:
			s1 = (s1-s0)//2
			x0 = xADD(x1, x0, xmin)
			x1 = xDBL(x1, ProjA)
		elif s1 % 2 == 0:
			s1 = s1//2
			xmin = xADD(x1, xmin, x0)
			x1 = xDBL(x1, ProjA)
		else:
			s0 = s0//2
			xmin = xADD(x0, xmin, x1)
			x0 = xDBL(x0, ProjA)

	while s1 % 2 == 0:
		s1 = s1//2
		x1 = xDBL(x1, ProjA)

	# Extra step with tripling as described in eprint 2017/212
	while s1 % 3 == 0:
		s1 = s1//3
		x1 = xTPL(x1, ProjA)

	if s1 > 1:
		x1 = xMUL(x1, s1, ProjA)

	return x1

#####################
# THREE-POINT LADDER
#####################

def ladder_3pt(P, Q, PmQ, m, ProjA):
	"""
	Computes 3pt-ladder (from SIKE)
	
	Input: projective points P, Q, PmQ = P-Q
		   scalar m, curve constant A = ProjA[0]/ProjA[1]
	Output: P+[m]Q
	"""

	binary = lambda n: n>0 and [n&1]+binary(n>>1) or []
	bits = binary(m)

	P0 = [Q[0], Q[1]]
	P1 = [P[0], P[1]]
	P2 = [PmQ[0], PmQ[1]]

	for i in range(len(bits)):
		if bits[i] == 1:
			P0, P1 = xDBLADD(P0, P1, P2, ProjA)
		else:
			P0, P2 = xDBLADD(P0, P2, P1, ProjA)

	return P1

def ladder_3ptaffs(P, Q, PmQ, m, ProjA):
	"""
	Computes 3pt-ladder (from SIKE)
	
	Input: projective points P, Q, PmQ = P-Q
		   scalar m, curve constant A = ProjA[0]/1
	Output: P+[m]Q
	"""

	assert ProjA[1] == [1,0]

	binary = lambda n: n>0 and [n&1]+binary(n>>1) or []
	bits = binary(m)

	P0 = [Q[0], Q[1]]
	P1 = [P[0], P[1]]
	P2 = [PmQ[0], PmQ[1]]

	for i in range(len(bits)):
		if bits[i] == 1:
			P0, P1 = xDBLADDaffs(P0, P1, P2, ProjA)
		else:
			P0, P2 = xDBLADDaffs(P0, P2, P1, ProjA)

	return P1




