from sage.all import *
import fpylll
from fpylll import IntegerMatrix, CVP
from fpylll.fplll.gso import MatGSO


###############################################
#                                             #
# Enumerate close vectors, from LearningToSQI #
#                                             #
###############################################

def solve_closest_vector_problem(lattice_basis, target):
    """
    Use the fpylll library to solve the CVP problem for a given
    lattice basis and target vector
    """
    L = IntegerMatrix.from_matrix(lattice_basis.LLL())
    v = CVP.closest_vector(L, target)
    # fpylll returns a type `tuple` object
    return vector(v)


def generate_short_vectors_fpyll(L, bound, count=2000):
    """
    Helper function for GenerateShortVectors and 
    generate_small_norm_quat which builds an iterator
    for short norm vectors of an LLL reduced lattice
    basis.
    """
    # # Move from Sage world to Fypll world
    A = IntegerMatrix.from_matrix(L)

    # Gram-Schmidt Othogonalization
    G = MatGSO(A)
    _ = G.update_gso()

    # Enumeration class on G with `count`` solutions
    # BEST_N_SOLUTIONS:
    # Starting with the nr_solutions-th solution, every time a new solution is found
    # the enumeration bound is updated to the length of the longest solution. If more
    # than nr_solutions were found, the longest is dropped.
    E = fpylll.Enumeration(
        G, nr_solutions=count, strategy=fpylll.EvaluatorStrategy.BEST_N_SOLUTIONS
    )

    # We need the row count when we call enumerate
    r = L.nrows()

    # If enumerate finds no solutions it raises an error, so we
    # wrap it in a try block
    try:
        # The arguments of enumerate are:
        # E.enumerate(first_row, last_row, max_dist, max_dist_expo)
        short_vectors = E.enumerate(0, r, bound, 0)
    except Exception as e:
        short_vectors = []
        
    return short_vectors

def generate_short_vectors(lattice_basis, bound, count=2000):
    """
    Generate a generator of short vectors with norm <= `bound`
    returns at most `count` vectors.
    
    Most of the heavy lifting of this function is done by 
    generate_short_vectors_fpyll
    """
    L = lattice_basis.LLL()
    short_vectors = generate_short_vectors_fpyll(L, bound, count=count)
    for _, xis in short_vectors:
        # Returns values x1,x2,...xr such that
        # x0*row[0] + ... + xr*row[r] = short vector
        v3 = vector([ZZ(xi) for xi in xis])
        v = v3 * L
        yield v


def generate_close_vectors(lattice_basis, target, p, L, count=2000):
    """
    Generate a generator of vectors which are close, without
    bound determined by N to the `target`. The first
    element of the list is the solution of the CVP.
    """
    # Compute the closest element
    closest = solve_closest_vector_problem(lattice_basis, target)
    yield closest

    # Now use short vectors below a bound to find
    # close enough vectors

    # Set the distance
    diff = target - closest
    distance = diff.dot_product(diff)

    # Compute the bound from L
    b0 = L // p
    bound = floor((b0 + distance) + (2 * (b0 * distance).sqrt()))

    short_vectors = generate_short_vectors(lattice_basis, bound, count=count)

    for v in short_vectors:
        yield closest + v


###############################################
#                                             #
#        Some basic helper functions          #
#                                             #
###############################################

def QuaternionOrderBasis(alpha, O):
    assert alpha in O
    Obasis = O.basis()
    M_O = Matrix(QQ, [b.coefficient_tuple() for b in Obasis]).transpose()
    vec_alpha = vector(alpha.coefficient_tuple())
    coeffs = M_O.inverse() * vec_alpha
    coeffs = [ZZ(c) for c in coeffs]
    #assert alpha == sum([c * beta for c, beta in zip(coeffs, Obasis)])
    return coeffs

def IdealGenerator(I):
    p = I.quaternion_algebra().ramified_primes()[0]
    bound = max(ceil(10*log(p,2)), 100)
    while True:
        alpha = sum([b * randint(-bound, bound) for b in I.basis()])
        if gcd(alpha.reduced_norm(), I.norm()**2) == I.norm():
            return alpha
        
def InverseIdeal(I):
    return I.conjugate()*(1/I.norm())

def MakePrimitive(alpha, O):
    v_alpha = QuaternionOrderBasis(alpha, O)
    d = gcd(v_alpha)
    assert alpha/d in O
    return alpha/d

def MakeCyclic(I):
    O = I.left_order()
    d = gcd([gcd(QuaternionOrderBasis(beta, O)) for beta in I.basis()])
    return I*(1/d)

def ReducedBasis(I):
    """
    Computes the Minkowski reduced basis of the
    input ideal. From LearningToSQI
    """
    def _matrix_to_gens(M, B):
        """
        Converts from a matrix to generators in the quat.
        algebra
        """
        return [sum(c * g for c, g in zip(row, B)) for row in M]

    B = I.basis()
    G = I.gram_matrix()
    U = G.LLL_gram().transpose()
    reduced_basis_elements = _matrix_to_gens(U, B)
    return reduced_basis_elements

def RandomEquivalentPrimeIdeal(I):
    p = I.quaternion_algebra().ramified_primes()[0]
    reduced_basis = ReducedBasis(I)
    bound = p**(0.5)*100 #100 is a bit arbitrary here...
    for _ in range(1000):
        coeffs = [randint(-10, 10) for _ in range(4)]
        delta = sum([c*beta for c, beta in zip(coeffs, reduced_basis)])
        dnorm = ZZ(delta.reduced_norm() / I.norm())
        if is_pseudoprime(dnorm) and dnorm < bound:
            return I*(delta.conjugate()/I.norm())
    return None

def IdealEquivalence(I, J):
    assert I.left_order() == J.left_order()
    IJ = I.conjugate() * J
    alpha = ReducedBasis(IJ)[0]/I.norm()
    assert J == I*alpha
    return alpha

def ConnectingIdeal(O0, O1):
    I = O0*O1
    I = I*I.norm().denominator()
    return I

###############################################
#                                             #
#           Basic ideal operations            #
#                                             #
###############################################

def pushforward(J, I):
    assert J.left_order() == I.left_order()
    return InverseIdeal(J)*J.intersection(I)

def pullback(J, I):
    assert J.right_order() == I.left_order()
    return J*I + J.left_order()*I.norm()


###############################################
#                                             #
#    Jumping between quaternion algebras      #
#                                             #
###############################################

def IsomorphismGamma(B_old, B):
    r"""
    Used for computing the isomorphism between B_old and B
    See Lemma 10 [EPSV23]
    """
    if B_old == B:
        return 1
    i_old, j_old, k_old = B_old.gens()
    q_old = -ZZ(i_old**2)
    i, j, k = B.gens()
    q = -ZZ(i**2) 
    p = -ZZ(j**2)
    x, y = DiagonalQuadraticForm(QQ, [1,p]).solve(q_old/q)
    return x + j*y, (x + j_old*y)**(-1)

def EvalIsomorphism(alpha, B, gamma):
    r"""
    Given alpha \in B_old, and gamma deteremining the isomorphism from B_old to B,
    returns alpha \in B
    """
    i, j, k = B.gens()
    return sum([coeff*b for coeff, b in zip(alpha.coefficient_tuple(), [1, i*gamma, j, k*gamma])]) 