from ec_misc import *
from fp_arith import get_p
from random import randint

####################################################################################################################
##################################          This file contains functions to       ##################################
##################################             compute the torsion basis          ##################################
####################################################################################################################


################################################
# SIKE techniques
################################################

def get_alpha(A):
	"""
	Function to find a root, alpha, of x^2 + Ax + 1 (Montgomery form of the Elliptic Curve)
	
	Input: A
	Output: alpha, a root of x^2 + Ax + 1
	"""
	assert A[1] == [1,0]
	inv2 = fp2_inv([2,0])

	alpha = fp2_sqr(A[0])
	alpha = fp2_sub(alpha, [4,0])
	alpha = fp2_squareroot(alpha)
	alpha = fp2_add(alpha, fp2_neg(A[0]))
	alpha = fp2_mul(alpha, [inv2, 0])       #can be better
	return alpha

def has_full_two_over_fp(Q, A, alpha):
	"""
	Function to determine if point Q has full 2-torsion on y^2 = x(x^2 + Ax + 1) 
	
	Input: Point Q, A, alpha s.t. x^2 + Ax + 1 = (x - alpha)(x - 1/alpha)
	Output: True,  if point Q has full 2-torsion and False otherwise
	"""
	assert Q[1] == [1,0]
	assert A[1] == [1,0]

	# Over Fp, Q[0] is always a square so we check if Q[0] - alpha is a non-square
	return not fp2_issquare(fp2_sub(Q[0], alpha))


################################################
# Torsion basis computation
################################################

def basis_two_torsion(ProjA, e):
	"""
	Function that outputs a basis of the 2^e-torsion on curve with constant A
 
	Input: affine curve constant ProjA = [A, [1,0]]
	  		power of two 2^e
	Output: points P, Q that form a basis of the 2^e-torsion
			difference PmQ = P-Q
	"""

	assert(ProjA[1] == [1,0])

	p = get_p()

	x = [1, randint(0,2**e)]		#to avoid peaks bug

	#find first point
	while True:		
		x = fp2_add(x, [0, 1])

		if not fp2_issquare(montgomery_rhs(x, ProjA[0])):
			continue

		#remove odd torsion
		m = (p+1)//(2**e)
		P = xMUL([x, [1, 0]], m, ProjA)

		#check for 2-torsion
		P2 = xDBLe(P, ProjA, e-1)
		if P2[1] != [0, 0]:
			break
		else:
			continue
	#find second point
	while True:
		x = fp2_add(x, [0, 1])

		if not fp2_issquare(montgomery_rhs(x, ProjA[0])):
			continue

		#remove odd torsion
		m = (p+1)//(2**e)
		Q = xMUL([x, [1, 0]], m, ProjA)

		#check for 2-torsion
		Q2 = xDBLe(Q, ProjA, e-1)

		if Q2[1] == [0, 0]:
			continue

		if not proj_point_equal(P2, Q2):
			break
		else:
			continue

	P, Q = double_proj_normalize(P, Q)

	PmQ = point_difference(P, Q, ProjA)

	return P, Q, PmQ


def basis_two_torsion_SIKE(ProjA, e):
	"""
	Function that outputs a basis of the 2^e-torsion on curve with constant A
	Follows SIKE spec
 
	Input: affine curve constant ProjA = [A, [1,0]]
	  		power of two 2^e
	Output: points P, Q that form a basis of the 2^e-torsion
			difference PmQ = P-Q
	"""	
 
	assert(ProjA[1] == [1,0])
	p = get_p()

	i = [1,0]
	x = i

	while True:
		x = fp2_add(x, i)

		if not fp2_issquare(montgomery_rhs(x, ProjA[0])):
			continue

		#remove odd torsion		
		m = (p+1)//(2**e)
		P = xMUL([x, [1, 0]], m, ProjA)

		#check for 2-torsion
		P2 = xDBLe(P, ProjA, e-1)
		if P2[1] != [0, 0]:
			break
		else:
			continue

	#find second point
	i = 1
	while True:
		i += 1
	
		x = random_nonsquare(i)

		if not fp2_issquare(montgomery_rhs(x, ProjA[0])):
			continue

		#remove odd torsion
		m = (p+1)//(2**e)
		Q = xMUL([x, [1, 0]], m, ProjA)

		#check for 2-torsion
		Q2 = xDBLe(Q, ProjA, e-1)

		if Q2[1] == [0, 0]:
			continue

		if not proj_point_equal(P2, Q2):
			break
		else:
			continue

	P, Q = double_proj_normalize(P, Q)

	PmQ = point_difference(P, Q, ProjA)

	return P, Q, PmQ

def basis_lam_torsion_APRES(ProjA, lam, f2):
	"""
	Function that outputs a basis of the 2^e-torsion on curve with constant A
	As done for ApresSQI

	Input: affine curve constant ProjA = [A, [1,0]]
			power of two 2^e
	Output: points P, Q that form a basis of the 2^e-torsion
			difference PmQ = P-Q
	"""	

	assert(ProjA[1] == [1,0])
	p = get_p()

	i = [1,0]
	x = i

	while True:
		x = fp2_add(x, i)

		if not fp2_issquare(montgomery_rhs(x, ProjA[0])):
			continue

		#remove odd torsion		
		m = (p+1)//(2**f2)
		P = xMULc([x, [1, 0]], m, ProjA)

		#then check amount of torsion p has, should always be enough (but unfortunate we have to check, try tricks)
		tor_p = 1
		P2 = xDBLaff(P, ProjA)
		prev_P = [P, P2]
		while P2[1] != [0,0]:
			tor_p += 1
			P2 = xDBLaff(P2, ProjA)
			prev_P.append(P2)

		if tor_p < lam:		#this almost never happens
			continue

		P = prev_P[tor_p - lam]			#this can be more elegant
		break

	#find second point
	i = 1
	while True:
		i += 1
		x = small_nonsquare(i)

		if not fp2_issquare(montgomery_rhs(x, ProjA[0])):
			continue

		Q = [x, [1,0]]
		#remove odd torsion
		# m = (p+1)//(2**f2)
		# Q = xMULc([x, [1, 0]], m, ProjA)
		# Q = xDBLe(Q, ProjA, f2 - lam)

		break
		# this is implicit
		# #check for 2-torsion
		# Q2 = xDBLe(Q, ProjA, e-1)

		# if Q2[1] == [0, 0]:
		# 	continue

		# if not proj_point_equal(P2, Q2):
		# 	break
		# else:
		# 	continue

	#P, Q = double_proj_normalize(P, Q)
	# PmQ = point_difference(P, Q, ProjA)

	return P, Q  #, PmQ


def basis_two_torsion_seed(ProjA, e):

	"""
	Function that outputs a basis of the 2^e-torsion on curve with constant A
	With seeding
 
	Input: affine curve constant ProjA = [A, [1,0]]
	  		power of two 2^e
	Output: points P, Q that form a basis of the 2^e-torsion,
			difference PmQ = P-Q, and seeds for P and Q
	"""

	assert ProjA[1] == [1,0]
	p = get_p()
	i = 1
	x = small_nonsquare(i)

	#find first point
	while True:
		i += 1
		x = small_nonsquare(i)

		if not fp2_issquare(montgomery_rhs(x, ProjA[0])):
			continue

		#remove odd torsion
		m = (p+1)//(2**e)
		P = xMUL([x, [1, 0]], m, ProjA)

		#check for 2-torsion
		P2 = xDBLe(P, ProjA, e-1)
		if P2[1] != [0, 0]:
			u0 = i
			break
		else:
			continue
	#find second point
	while True:
		i += 1
		x = small_nonsquare(i)

		if not fp2_issquare(montgomery_rhs(x, ProjA[0])):
			continue

		#remove odd torsion
		m = (p+1)//(2**e)
		Q = xMUL([x, [1, 0]], m, ProjA)

		#check for 2-torsion
		Q2 = xDBLe(Q, ProjA, e-1)

		if Q2[1] == [0, 0]:
			continue

		if not proj_point_equal(P2, Q2):
			u1 = i
			break
		else:
			continue

	P, Q = double_proj_normalize(P, Q)

	PmQ = point_difference(P, Q, ProjA)

	return P, Q, PmQ, [u0,u1]

def basis_two_torsion_APRES(ProjA, e, u):

	"""
	Function that outputs a basis of the 2^e-torsion on curve with constant A
	As done for ApresSQI
 
	Input: affine curve constant ProjA = [A, [1,0]]
	  		power of two 2^e
	Output: points P, Q that form a basis of the 2^e-torsion
			difference PmQ = P-Q
	"""
	assert(ProjA[1] == [1,0])
	p = get_p()

	x = small_nonsquare(u[0])

	#remove odd torsion
	m = (p+1)//(2**e)
	P = xMULc([x, [1, 0]], m, ProjA)

	x = small_nonsquare(u[1])
	Q = xMULc([x, [1, 0]], m, ProjA)

	P, Q = double_proj_normalize(P, Q)
	PmQ = point_difference(P, Q, ProjA)

	return P, Q, PmQ


def complete_basis_two_torsion(P, ProjA, e):

	"""
	Function that completes the basis for 2^e-torsion

	Input: affine point P of order 2^e
			affine curve constant ProjA = [A, [1,0]]
			power of two 2^e
	Output: point Q such that (P,Q) form a basis of the 2^e-torsion
			difference PmQ = P-Q
	"""

	assert(ProjA[1] == [1,0])
	assert(P[1] == [1,0])
	p = get_p()

	#compute point of order 2 from P
	P2 = xDBLe(P, ProjA, e-1)

	#find second point
	x = [1, randint(0,2**e)]		#to avoid peaks bug

	while True:
		x = fp2_add(x, [0, 1])

		if not fp2_issquare(montgomery_rhs(x, ProjA[0])):
			continue

		#remove odd torsion
		m = (p+1)//(2**e)
		Q = xMUL([x, [1, 0]], m, ProjA)

		#check for 2-torsion
		Q2 = xDBLe(Q, ProjA, e-1)

		if Q2[1] == [0, 0]:
			continue

		if not proj_point_equal(P2, Q2):
			break
		else:
			continue

	Q = proj_point_normalize(Q)
	PmQ = point_difference(P, Q, ProjA)

	return P, Q, PmQ

def complete_basis_two_torsion_SIKE(P, ProjA, e):

	"""
	Function that completes the basis for 2^e-torsion
	Finding the second point uses SIKE's algorithm a la 2016/963
	
	Input: affine point P of order 2^e
			affine curve constant ProjA = [A, [1,0]]
			power of two 2^e
	Output: point Q such that (P,Q) form a basis of the 2^e-torsion
			difference PmQ = P-Q
	"""

	assert(ProjA[1] == [1,0])
	assert(P[1] == [1,0])
	p = get_p()

	#find second point
	i = 1
	x = random_nonsquare(i)

	while True:
		i += 1
		x = random_nonsquare(i)

		if not fp2_issquare(montgomery_rhs(x, ProjA[0])):
			continue

		#remove odd torsion
		m = (p+1)//(2**e)
		Q = xMUL([x, [1, 0]], m, ProjA)

		break

	Q = proj_point_normalize(Q)
	PmQ = point_difference(P, Q, ProjA)

	return P, Q, PmQ

def complete_basis_two_torsion_seed_get(P, ProjA, e):
	"""
	Function that completes the basis for 2^e-torsion with seeding
	Finding the second point uses SIKE's algorithm a la 2016/963

	Input: affine point P of order 2^e
	 		affine curve constant ProjA = [A, [1,0]]
	 		power of two 2^e
	Output: point Q such that (P,Q) form a basis of the 2^e-torsion
			difference PmQ = P-Q and seed i that gives Q = [[i,0], [1,0]]
	"""

	assert(ProjA[1] == [1,0])
	assert(P[1] == [1,0])
	p = get_p()

	#find second point
	i = -1
	x = random_nonsquare(i)

	while True:
		i += 1
		x = random_nonsquare(i)

		if not fp2_issquare(montgomery_rhs(x, ProjA[0])):
			continue

		#remove odd torsion
		m = (p+1)//(2**e)
		Q = xMUL([x, [1, 0]], m, ProjA)

		break

	Q = proj_point_normalize(Q)
	PmQ = point_difference(P, Q, ProjA)

	return P, Q, PmQ, i


def complete_basis_two_torsion_seed_get_plus(ProjA, e):
	"""
	Function that completes the basis for 2^e-torsion with seeding
	Finding the second point uses SIKE's algorithm a la 2016/963

	Input: affine point P of order 2^e
	 		affine curve constant ProjA = [A, [1,0]]
	 		power of two 2^e
	Output: point Q such that (P,Q) form a basis of the 2^e-torsion
			difference PmQ = P-Q and 
			seeds i,j that give Q = [[i,0], [1,0]] and P = [[small_non_square(j)],[1,0]]

	"""
	
	assert(ProjA[1] == [1,0])
	p = get_p()
	#find first point
	i = -1
	while True:
		i += 1
		P = [small_nonsquare(i),[1,0]]

		if not fp2_issquare(montgomery_rhs(P[0], ProjA[0])):
			continue

		#remove odd torsion
		m = (p+1)//(2**e)
		P1 = xMUL(P, m, ProjA)

		#check for 2-torsion
		P2 = xDBLe(P1, ProjA, e-1)
		if P2[1] == [0, 0]:
			continue

		break

	#find second point
	j = -1

	while True:
		j += 1
		Q = [small_square_with_minus_one_nonsquare(j), [1,0]]

		if not fp2_issquare(montgomery_rhs(Q[0], ProjA[0])):
			continue
		
		#remove odd torsion
		m = (p+1)//(2**e)
		Q1 = xMUL(P, m, ProjA)

		#check for 2-torsion
		Q2 = xDBLe(Q1, ProjA, e-1)
		if Q2[1] == [0, 0]:
			continue

		break

	return P, Q, [i, j]


def complete_basis_chall_torsion_seed_get_plus(ProjA, e2, e3):
	"""
	Function that completes the basis for 2^e2*3^e3-torsion with seeding (used in Challenge phase)
	Finding the second point uses SIKE's algorithm a la 2016/963

	Input: affine point P of order 2^e2*3^e3
	 		affine curve constant ProjA = [A, [1,0]]
	 		power of two 2^e2 and 3^e3
	Output: point Q such that (P,Q) form a basis of the 2^e2*3^e3-torsion
			difference PmQ = P-Q and 
			seeds i,j that give Q = [[i,0], [1,0]] and P = [[small_non_square(j)],[1,0]]
	"""

	assert(ProjA[1] == [1,0])
	alpha = get_alpha(ProjA)
	p = get_p()
	m = (p+1)//(2**e2*3**e3)

	#find first point
	j = -1

	while True:
		j += 1
		Q = [small_square_with_minus_one_nonsquare(j),[1,0]]

		if not fp2_issquare(montgomery_rhs(Q[0], ProjA[0])):
			continue

		if not has_full_two_over_fp(Q, ProjA, alpha):
			continue

		#remove other torsion
		Q1 = xMUL(Q, m, ProjA)


		#check for 2 and 3-torsion
		Q23 = xTPLe(xDBLe(Q1, ProjA, e2-1), ProjA, e3-1)

		Q2 = xTPL(Q23, ProjA)
		if Q2[1] == [0,0]:
			continue
		Q3 = xDBL(Q23, ProjA)
		if Q3[1] == [0, 0]:
			continue
		# Point is valid
		break

	#find second point
	i = -1

	while True:
		i += 1
		P = [small_nonsquare(i), [1,0]]

		if not fp2_issquare(montgomery_rhs(P[0], ProjA[0])):
			continue

				#remove other torsion
		P1 = xMUL(P, m, ProjA)

		#check for 2 and 3-torsion
		P23 = xTPLe(xDBLe(P1, ProjA, e2-1), ProjA, e3-1)

		P2 = xTPL(P23, ProjA)
		if P2[1] == [0,0]:
			continue
		P3 = xDBL(P23, ProjA)
		if P3[1] == [0, 0]:
			continue
		
		#Check 2 basis
		if proj_point_equal(P2, Q2):
			continue

		#Check 3 basis
		if proj_point_equal(P3, Q3):
			continue
		if proj_point_equal(xDBL(P3, ProjA), Q3):
			continue

		#basis is valid

		break

	return P, Q, [i, j]

def complete_basis_two_torsion_seed(P, ProjA, e, u):
	"""
	Function that completes the basis for 2^e-torsion with seed for second point
	Finding the second point uses SIKE's algorithm a la 2016/963

	Input: affine point P of order 2^e
	 		affine curve constant ProjA = [A, [1,0]]
	 		power of two 2^e
			seed u
	Output: point Q such that (P,Q) form a basis of the 2^e-torsion
			difference PmQ = P-Q and 
	"""
	assert(ProjA[1] == [1,0])
	assert(P[1] == [1,0])
	p = get_p()

	#second point is given by seed
	x = random_nonsquare(u)

	#remove odd torsion
	m = (p+1)//(2**e)
	Q = xMUL([x, [1, 0]], m, ProjA)
	Q = proj_point_normalize(Q)
	PmQ = point_difference(P, Q, ProjA)

	return P, Q, PmQ

def basis_chall_torsion_apres(ProjA, e2, e3):
	"""
	Function that completes the basis for 2^e2*3^e3-torsion with seeding (used in Challenge phase)
	As done in ApresSQI

	Input: affine point P of order 2^e2*3^e3
	 		affine curve constant ProjA = [A, [1,0]]
	 		power of two 2^e2 and 3^e3
	Output: point Q such that (P,Q) form a basis of the 2^e2*3^e3-torsion
			difference PmQ = P-Q and 
			seeds i,j that give Q = [[i,0], [1,0]] and P = [[small_non_square(j)],[1,0]]

	"""    
	assert(ProjA[1] == [1,0])
	alpha = get_alpha(ProjA)
	p = get_p()
	m = (p+1)//(2**e2*3**e3)

	#find first point
	i = 1

	while True:
		i += 1
		Q = [[i,0],[1,0]]

		if not fp2_issquare(montgomery_rhs(Q[0], ProjA[0])):
			continue

		if not has_full_two_over_fp(Q, ProjA, alpha):
			continue

		#remove other torsion
		Q1 = xMUL(Q, m, ProjA)


		#check for 2 and 3-torsion
		Q23 = xTPLe(xDBLe(Q1, ProjA, e2-1), ProjA, e3-1)

		Q2 = xTPL(Q23, ProjA)
		if Q2[1] == [0,0]:
			continue
		Q3 = xDBL(Q23, ProjA)
		if Q3[1] == [0, 0]:
			continue
		# Point is valid
		break

	#find second point
	j = 1

	while True:
		j += 1
		P = [small_nonsquare(j), [1,0]]

		if not fp2_issquare(montgomery_rhs(P[0], ProjA[0])):
			continue

				#remove other torsion
		P = xMUL(P, m, ProjA)

		#check for 2 and 3-torsion
		P23 = xTPLe(xDBLe(P, ProjA, e2-1), ProjA, e3-1)

		P2 = xTPL(P23, ProjA)
		if P2[1] == [0,0]:
			continue
		P3 = xDBL(P23, ProjA)
		if P3[1] == [0, 0]:
			continue
		
		#Check 2 basis
		if proj_point_equal(P2, Q2):
			continue

		#Check 3 basis
		if proj_point_equal(P3, Q3):
			continue
		if proj_point_equal(xDBL(P3, ProjA), Q3):
			continue

		#basis is valid
		break

	return P, Q1



def basis_chal_torsion(ProjA, e2, e3):
	"""
	Function that completes the basis for 2^e2*3^e3-torsion (used in Challenge phase)
	Uses basis_two_torsion if e3 == 0 to find basis for 2^e2-torsion
 
	Input: affine curve constant ProjA = [A, [1,0]]
	 		power of two 2^e2 
			power of three 3^e3
	Output: point P,Q such that (P,Q) form a basis of the 2^e2*3^e3-torsion
			difference PmQ = P-Q and 

	"""
	if e3 == 0:
		#the improved version is available as basis_lam_torsion_APRES
		return basis_two_torsion(ProjA, e2) 

	assert(ProjA[1] == [1,0])
	p = get_p()

	x = [1, randint(0,2**e2)]		#to avoid peaks bug

	# to remove cofactors
	m = (p+1)//(2**e2*3**e3)
	#find first point
	while True:
		x = fp2_add(x, [0, 1])

		if not fp2_issquare(montgomery_rhs(x, ProjA[0])):
			continue

		#remove other torsion
		P = xMUL([x, [1, 0]], m, ProjA)


		#check for 2 and 3-torsion
		P23 = xTPLe(xDBLe(P, ProjA, e2-1), ProjA, e3-1)

		P2 = xTPL(P23, ProjA)
		if P2[1] == [0,0]:
			continue
		P3 = xDBL(P23, ProjA)
		if P3[1] == [0, 0]:
			continue
		# Point is valid
		break

	#find second point
	while True:
		x = fp2_add(x, [0, 1])

		if not fp2_issquare(montgomery_rhs(x, ProjA[0])):
			continue

		#remove other torsion
		Q = xMUL([x, [1, 0]], m, ProjA)

		#check for 2 and 3-torsion
		Q23 = xTPLe(xDBLe(Q, ProjA, e2-1), ProjA, e3-1)

		Q2 = xTPL(Q23, ProjA)
		if Q2[1] == [0,0]:
			continue
		Q3 = xDBL(Q23, ProjA)
		if Q3[1] == [0, 0]:
			continue

		#Check 2 basis
		if proj_point_equal(P2, Q2):
			continue

		#Check 3 basis
		if proj_point_equal(P3, Q3):
			continue
		if proj_point_equal(xDBL(P3, ProjA), Q3):
			continue

		#basis is valid
		break

	P, Q = double_proj_normalize(P, Q)

	PmQ = point_difference(P, Q, ProjA)

	return P, Q, PmQ


def basis_chal_torsion_SIKE(ProjA, e2, e3):

	"""
	Function that completes the basis for 2^e2*3^e3-torsion (used in Challenge phase)
	Uses basis_two_torsion_SIKE if e3 == 0 to find basis for 2^e2-torsion
 
	Input: affine curve constant ProjA = [A, [1,0]]
	 		power of two 2^e2 
			power of three 3^e3
	Output: point P,Q such that (P,Q) form a basis of the 2^e2*3^e3-torsion
			difference PmQ = P-Q and 

	"""

	if e3 == 0:
		#the improved version is available as basis_lam_torsion_APRES
		return basis_two_torsion_SIKE(ProjA, e2) 

	assert(ProjA[1] == [1,0])
	p = get_p()

	# to remove cofactors
	m = (p+1)//(2**e2*3**e3)
	#find first point
	i = 1
	x = random_nonsquare(i)

	while True:
		i += 1
		x = random_nonsquare(i)

		if not fp2_issquare(montgomery_rhs(x, ProjA[0])):
			continue

		#remove other torsion
		P = xMUL([x, [1, 0]], m, ProjA)


		#check for 2 and 3-torsion
		P23 = xTPLe(xDBLe(P, ProjA, e2-1), ProjA, e3-1)

		P2 = xTPL(P23, ProjA)
		if P2[1] == [0,0]:
			continue
		P3 = xDBL(P23, ProjA)
		if P3[1] == [0, 0]:
			continue
		# Point is valid
		break

	while True:
		i += 1
		x = random_nonsquare(i)

		if not fp2_issquare(montgomery_rhs(x, ProjA[0])):
			continue

		#remove other torsion
		Q = xMUL([x, [1, 0]], m, ProjA)

		#check for 2 and 3-torsion
		Q23 = xTPLe(xDBLe(Q, ProjA, e2-1), ProjA, e3-1)

		Q2 = xTPL(Q23, ProjA)
		if Q2[1] == [0,0]:
			continue
		Q3 = xDBL(Q23, ProjA)
		if Q3[1] == [0, 0]:
			continue

		#Check 2 basis
		if proj_point_equal(P2, Q2):
			continue

		#Check 3 basis
		if proj_point_equal(P3, Q3):
			continue
		if proj_point_equal(xDBL(P3, ProjA), Q3):
			continue

		#basis is valid
		break

	P, Q = double_proj_normalize(P, Q)

	PmQ = point_difference(P, Q, ProjA)

	return P, Q, PmQ


def basis_three_torsion(ProjA, e):
	"""
	Function that finds a basis P,Q for the 3^e-torsion

	Input: affine curve constant ProjA = [A, [1,0]]
			power of three 3^e
	Output: points P, Q that form a basis of the 3^e-torsion 
			difference PmQ = P-Q
	"""

	assert(ProjA[1] == [1,0])
	p = get_p()

	x = [1, randint(0,2**e)]		#to avoid peaks bug

	#find first point
	while True:
		x = fp2_add(x, [0, 1])

		if not fp2_issquare(montgomery_rhs(x, ProjA[0])):
			continue

		#remove other torsion
		m = (p+1)//(3**e)
		P = xMUL([x, [1, 0]], m, ProjA)

		#check for 3-torsion
		P3 = xTPLe(P, ProjA, e-1)

		if P3[1] != [0, 0]:
			break
		else:
			continue

	#find second point
	while True:
		x = fp2_add(x, [0, 1])

		if not fp2_issquare(montgomery_rhs(x, ProjA[0])):
			continue

		#remove other torsion
		m = (p+1)//(3**e)
		Q = xMUL([x, [1, 0]], m, ProjA)

		#check for 3-torsion
		Q3 = xTPLe(Q, ProjA, e-1)

		if Q3[1] == [0, 0]:
			continue

		if not proj_point_equal(P3, Q3):
			break
		else:
			continue

	P, Q = double_proj_normalize(P, Q)

	PmQ = point_difference(P, Q, ProjA)

	return P, Q, PmQ


def complete_basis_three_torsion(P, ProjA, e):
	"""
	Function that completes the basis for the 3^e-torsion (given P, finds another point Q)

	Input: affine point P of order 3^e
 			affine curve constant ProjA = [A, [1,0]]
			power of three 3^e
	Output: points P, Q that form a basis of the 3^e-torsion 
			difference PmQ = P-Q
	"""

	assert(ProjA[1] == [1,0])
	assert(P[1] == [1,0])
	p = get_p()

	#compute point of order 3 from P
	P3 = xTPLe(P, ProjA, e-1)

	#find second point
	x = [1, randint(0,2**e)]		#to avoid peaks bug

	while True:
		x = fp2_add(x, [0, 1])

		if not fp2_issquare(montgomery_rhs(x, ProjA[0])):
			continue

		#remove other torsion
		m = (p+1)//(3**e)
		Q = xMUL([x, [1, 0]], m, [ProjA[0], [1, 0]])

		#check for 3-torsion
		Q3 = xTPLe(Q, [ProjA[0], [1, 0]], e-1)

		if Q3[1] == [0, 0]:
			continue

		if not proj_point_equal(P3, Q3):
			break
		else:
			continue

	Q = proj_point_normalize(Q)
	PmQ = point_difference(P, Q, ProjA)

	return P, Q, PmQ


def point_difference(P, Q, ProjA):
	"""
	Function that computes the difference P-Q given two points P,Q lying on curve with projective constant ProjA
	Follows code from SQIsign NIST version 1.0

	Input: Affine points P,Q
			Affine curve constant ProjA = [A, [1,0]]
	Output: x-coordinate of P-Q = PmQ
 
	"""

	#check if all inputs are affine
	if ProjA[1] != [1,0] or P[1] != [1,0] or Q[1] != [1,0]:
		print_error("Input to point_difference must be affine!")

	PmQZ = fp2_sub(P[0], Q[0])
	t2 = fp2_mul(P[0], Q[0])
	t3 = fp2_sub(t2, [1,0])
	t0 = fp2_mul(PmQZ, t3)
	PmQZ = fp2_sqr(PmQZ)
	t0 = fp2_sqr(t0)
	t1 = fp2_add(t2, [1,0])
	t3 = fp2_add(P[0], Q[0])
	t1 = fp2_mul(t1, t3)
	t2 = fp2_mul(t2, ProjA[0])
	t2 = fp2_add(t2, t2)
	t1 = fp2_add(t1, t2)
	t2 = fp2_sqr(t1)
	t0 = fp2_sub(t2, t0)
	assert fp2_issquare(t0)		#TODO: is this required or can we save it?
	t0 = fp2_squareroot(t0)
	PmQX = fp2_add(t0, t1)

	return [PmQX, PmQZ]



##################################
#
#      Improved, see Thm. 3
#
##################################

"""
def sample_initial_Q(A, m, f2):
    #given affine A mont coeff
    #samples points from Fp until it has 2^f torsion
    #this implies it must be over (0,0)
    #returns the point with order 2^f

    #should probably not be used!!! instead use improved sampling below

    assert A[1] == [1,0]

    i = 1

    while True:
        i += 1
        Q = [[i,0],[1,0]]

        if not fp2_issquare(montgomery_rhs(Q[0], A[0])):
            continue

        #remove odd torsion
        Q = xMULc(Q, m, A)

        #check for 2-torsion
        Q2 = xDBLe_aff(Q, A, f2 - 1)
        if Q2[1] == [0, 0]:
            continue

        assert Q2[0] == [0,0]

        break

    return Q"""

def sample_improved_Q(A):
    """
		Given affine Montgomery coefficient A, samples points from Fp until it has order 2^f
		This implies it must be over (0,0)
		Returns a point with order 2^f
	"""

    assert A[1] == [1,0]
    alpha = get_alpha(A)

    i = 1

    while True:
        i += 1
        Q = [[i,0],[1,0]]

        if not fp2_issquare(montgomery_rhs(Q[0], A[0])):
            continue

        if not has_full_two_over_fp(Q, A, alpha):
            continue

        break

    return Q

def sample_improved_Q_plus_plus(A):
    """
		As above but now we use the fact that alpha must be square
  		Hence, given n in Fp2 that is square such that n-1 is nonsquare
    	x = n*alpha will be square, but x - alpha = (n-1)*alpha is nonsquare
	"""

    assert A[1] == [1,0]
    alpha = get_alpha(A)

    i = 0

    while True:
        n = small_square_with_minus_one_nonsquare(i)
        Q = [fp2_mul(n, alpha), [1,0]]
        i += 1

        if not fp2_issquare(montgomery_rhs(Q[0], A[0])):
            continue

        #assert has_full_two_over_fp(Q, A, alpha)

        break

    return Q
    
def sample_secondary_P(A):
    j = 0
    while True:
        P = [small_nonsquare(j), [1,0]]
        j += 1

        if not fp2_issquare(montgomery_rhs(P[0], A[0])):
            continue

        break

    return P