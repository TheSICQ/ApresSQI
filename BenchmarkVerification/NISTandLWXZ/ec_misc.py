
####################################################################################################################
##################################          This file contains miscellaneous      ##################################
##################################             elliptic curve functions           ##################################
####################################################################################################################

from mont_xonly import *
from random import randint
from fp_arith import get_p

################################################
# EC-related functions/conversions
################################################

def j_invariant(ProjA):
	"""
	Function that computes the j-invariant of curve with Montgomery curve constant ProjA

	Input: curve contant A = [A,C] = ProjA[0]/ProjA[1]
	Output: j-invariant
	"""

	j = fp2_sqr(ProjA[0])
	t1 = fp2_sqr(ProjA[1])
	t0 = fp2_add(t1, t1)
	t0 = fp2_sub(j,t0)
	t0 = fp2_sub(t0,t1)
	j = fp2_sub(t0,t1)
	t1 = fp2_sqr(t1)
	j = fp2_mul(j, t1)
	t0 = fp2_add(t0, t0)
	t0 = fp2_add(t0, t0)
	t1 = fp2_sqr(t0)
	t0 = fp2_mul(t0, t1)
	t0 = fp2_add(t0, t0)
	t0 = fp2_add(t0, t0)
	j = fp2_inv(j)
	j = fp2_mul(t0, j)

	return j


def ProjA_to_ProjACplus(ProjA):
	"""
	Function that converts Montgomery curve constant ProjA = [A, C] to ProjA24Plus = [A+2C, 4C]

	Input: curve contant A = [A,C] = ProjA[0]/ProjA[1]
	Output: ProjAPlus = [A24plus, C24] = [A+2C, 4C]
	"""

	C24 = fp2_add(ProjA[1], ProjA[1])
	A24plus = fp2_add(ProjA[0], C24)
	C24 = fp2_add(C24, C24)

	return [A24plus, C24]


def ProjACplus_to_ProjA(ProjAPlus):
	"""
	Function that converts curve constant ProjA24Plus = [A+2C, 4C] to ProjA = [A, C]

	Input: ProjAPlus = [A24plus, C24] = [A+2C, 4C]
	Output: A = [4A,4C] = ProjA[0]/ProjA[1]
	"""

	A = fp2_add(ProjAPlus[0], ProjAPlus[0])
	A = fp2_sub(A, ProjAPlus[1])
	A = fp2_add(A, A)

	return [A, ProjAPlus[1]]


def ProjACplusminus_to_ProjA(A24plusminus):
	"""
	Function that converts curve constant ProjACplusminus = [A+2C, A-2C] to ProjA = [A, C]

	Input: ProjAPlus = [A24plus, A24minus] = [A+2C, A-2C]
	Output: A = [4A,4C] = ProjA[0]/ProjA[1]
	"""

	A = fp2_add(A24plusminus[0], A24plusminus[1]) 	# 2A
	A = fp2_add(A, A) 								# 4A
	C = fp2_sub(A24plusminus[0], A24plusminus[1]) 	# 4C

	return [A, C]


def ProjA_to_ProjACplusminus(ProjA):
	"""
	Function that converts curve constant ProjA = [A, C] to ProjACplusminus = [A+2C, A-2C]

	Input: ProjA = [A,C] = ProjA[0]/ProjA[1]
	Output: ProjAPlus = [A24plus, A24minus] = [A+2C, A-2C]
	"""

	Aminus = fp2_add(ProjA[1], ProjA[1]) 	# 2C
	Aplus = fp2_add(ProjA[0], Aminus) 		# A + 2C
	Aminus = fp2_sub(ProjA[0], Aminus) 		# A - 2C

	return [Aplus, Aminus]


def proj_point_equal(P, Q):
	"""
	Function that determines if two projective points represent the same affine point

	Input: two projective points P,Q
	Output: True if P, Q represent the same affine point, False otherwise
	"""

	if P == [[0,0],[0,0]] or Q == [[0,0],[0,0]]:
		return P == Q

	t0 = fp2_mul(P[0], Q[1])
	t1 = fp2_mul(P[1], Q[0])

	return t0 == t1


def proj_point_normalize(Pp):
	"""
	Function that normalises a projective point

	Input: projective point P
	Output: affine representation of P
	"""

	P = [fp2_copy(p) for p in Pp] # Pp.copy() gives errors 

	P[1] = fp2_inv(P[1])
	P[0] = fp2_mul(P[0], P[1])
	P[1] = [1,0]

	return P

def double_proj_normalize(Pp, Qq):
	"""
	Function that normalises two projective point

	Input: projective points P, Q
	Output: affine representations of P, Q
	"""
 
	P = [fp2_copy(p) for p in Pp] # Pp.copy() gives errors 
	Q = [fp2_copy(q) for q in Qq] # Qq.copy() gives errors 


	R = fp2_mul(P[1], Q[1])
	R = fp2_inv(R)
	P[0] = fp2_mul(P[0], Q[1])
	P[0] = fp2_mul(P[0], R)
	Q[0] = fp2_mul(Q[0], P[1])
	Q[0] = fp2_mul(Q[0], R)
	P[1] = [1,0]
	Q[1] = [1,0]

	return P, Q

def montgomery_rhs(x, A):
	"""
	Function that evaluates the curve equation RHS x^3+Ax^2+x for a given x
 
	Input: x and curve constant A
	Output: x^3+Ax^2+x
	"""
	rhs = fp2_sqr(x)
	t0 = fp2_mul(A, x)
	rhs = fp2_add(rhs, t0)  #x^2+Ax
	rhs = fp2_add(rhs, [1, 0])  #x^2+Ax+1
	rhs = fp2_mul(rhs, x)
	return rhs


def normalize_curve(oldProjA):
	"""
	Function that computes an isomorphism from a curve to a canonical representative of the isomorphism class
	Follows code from SQIsign NIST version 1.0
 
	Input: projective curve parameter oldProjA = (A : C)
	Output: new curve parameter, isomorphism
	"""

	t0 = fp2_sqr(oldProjA[1])
	t1 = fp2_add(t0, t0)
	t2 = fp2_add(t1, t1)
	t3 = fp2_sqr(oldProjA[0])
	t2 = fp2_sub(t3, t2)
	t2 = fp2_squareroot(t2)
	t0 = fp2_add(t0, t1)
	t1 = fp2_mul(t2, t1)
	t5 = fp2_sub(t3, t0)
	t5 = fp2_mul(t5, oldProjA[0])
	t4 = fp2_add(t0, t0)
	t0 = fp2_add(t4, t0)
	t0 = fp2_sub(t0, t3)
	t3 = fp2_add(t3, t3)
	t3 = fp2_mul(t3, t2)
	t2 = fp2_mul(t2, t0)
	t0 = fp2_add(t2, t5)
	t2 = fp2_sub(t2, t5)
	t1 = fp2_inv(t1)
	t0 = fp2_mul(t0, t1)
	t2 = fp2_mul(t2, t1)
	t1 = fp2_mul(t3, t1)

	#choose some canonical solution
	#here: smallest real part
	if t0[0] > t1[0] or (t0[0] == t1[0] and t0[1] > t1[1]):
		t0 = fp2_copy(t1)
	if t0[0] > t2[0] or (t0[0] == t2[0] and t0[1] > t2[1]):
		t0 = fp2_copy(t2)

	newA = fp2_squareroot(t0)
	newC = [1, 0]
	isom = ec_isomorphism(oldProjA, [newA, newC])

	return [newA, newC], isom


def ec_isomorphism(fromProjA, toProjA):
	"""
	Function that computes an isomorphism from a curve to another curve of the isomorphism class
	Follows code from SQIsign NIST version 1.0

	Input: domain curve fromProjA, codomain curve toProjA
	Output: isomorphism data isom = [isomD, isomNx, isomNz]
	"""

	t0 = fp2_mul(fromProjA[0], toProjA[1])
	t0 = fp2_sqr(t0)
	t1 = fp2_mul(toProjA[0], fromProjA[1])
	t1 = fp2_sqr(t1)
	t2 = fp2_mul(toProjA[1], fromProjA[1])
	t2 = fp2_sqr(t2)
	t3 = fp2_add(t2, t2)
	t2 = fp2_add(t3, t2)
	t3 = fp2_sub(t2, t0)
	t4 = fp2_sub(t2, t1)
	t3 = fp2_inv(t3)
	t4 = fp2_mul(t4, t3)
	t4 = fp2_squareroot(t4)
	t3 = fp2_sqr(t4)
	t3 = fp2_mul(t3, t4)

	t0 = fp2_sqr(fromProjA[1])
	t1 = fp2_add(t0, t0)
	t1 = fp2_add(t1, t1)
	t1 = fp2_add(t1, t1)
	t0 = fp2_add(t0, t1)
	t2 = fp2_sqr(fromProjA[0])
	t2 = fp2_add(t2, t2)
	t2 = fp2_sub(t2, t0)
	t2 = fp2_mul(t2, fromProjA[0])
	t0 = fp2_sqr(toProjA[1])
	t0 = fp2_mul(t0, toProjA[1])
	t2 = fp2_mul(t2, t0)
	t3 = fp2_mul(t3, t2)
	t0 = fp2_sqr(toProjA[1])
	t1 = fp2_add(t0, t0)
	t1 = fp2_add(t1, t1)
	t1 = fp2_add(t1, t1)
	t0 = fp2_add(t0, t1)
	t2 = fp2_sqr(toProjA[0])
	t2 = fp2_add(t2, t2)
	t2 = fp2_sub(t2, t0)
	t2 = fp2_mul(t2, toProjA[0])
	t0 = fp2_sqr(fromProjA[1])
	t0 = fp2_mul(t0, fromProjA[1])
	t2 = fp2_mul(t2, t0)
	if t2 != t3:
		t4 = fp2_neg(t4)

	t0 = [1, 0]
	isomD = fp2_add(t0, t0)
	isomD = fp2_add(isomD, t0)
	isomD = fp2_mul(isomD, fromProjA[1])
	isomD = fp2_mul(isomD, toProjA[1])
	isomNx = fp2_mul(isomD, t4)
	t4 = fp2_mul(t4, fromProjA[0])
	t4 = fp2_mul(t4, toProjA[1])
	t0 = fp2_mul(toProjA[0], fromProjA[1])
	isomNz = fp2_sub(t0, t4)

	return [isomD, isomNx, isomNz]


def ec_isomorphism_eval(P, isom):
	"""
	Function that evaluates an isomorphism at a point P
	Follows code from SQIsign NIST version 1.0

	Input: point P to evaluate, isomorphism data isom = [isomD, isomNx, isomNz]
	Output: evaluation of P = [EvalPX, EvalPZ]
	"""

	EvalPX = fp2_mul(P[0], isom[1])
	t0 = fp2_mul(P[1], isom[2])
	EvalPX = fp2_sub(EvalPX, t0)
	EvalPZ = fp2_mul(P[1], isom[0])

	return [EvalPX, EvalPZ]

def is_supersingular(ProjA):
	"""
	Determines whether a curve, defined by the constant ProjA, is supersingular 
 
	Input: projective montogery coefficient ProjA
	Output: True or False 
				If False: definitely ordinary
				If True: Supersingular with very high probability
	"""

	A = proj_point_normalize(ProjA)
	p = get_p()
	for i in range(10):
		x = [randint(0,p), 0]		#to avoid peaks bug
		if not fp2_issquare(montgomery_rhs(x, A[0])):
			continue

		P = xMUL([x, [1, 0]], p+1, A)

		if P[1] != [0,0]:
			return False
	
	return True