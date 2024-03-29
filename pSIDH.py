from itertools import chain, combinations
from sage.all import CRT, DLP, GF, ProjectiveSpace, ZZ, Zmod, choice, factor, gcd, inverse_mod, mod, prod

from deuring import randomideal
from sqisign.deuring import IdealToIsogenyFromKLPT
from sqisign.KLPT import EichlerModConstraint, EquivalentPrimeIdealHeuristic, IdealModConstraint, StrongApproximationHeuristic


# Finds a connecting ideal I between O1 and O2
# i.e., an O1-left ideal and O2-right ideal
def ConnectingIdeal(params, O1, O2):
    B = params['B']

    assert O1.discriminant() == B.discriminant(), 'O1 is not a maximal ideal'
    assert O2.discriminant() == B.discriminant(), 'O2 is not a maximal ideal'

    N = O1.intersection(O2).free_module().index_in(O1.free_module())  # Sage orders in quaternion algebras don't directly support index_in()
    I = N * O1 * O2

    assert I.left_order() == O1 and I.right_order() == O2, 'I is not an O1,O2-connecting ideal'

    return I


def EichlerSuborderNormEquation_helper(params, D, I, N):
    ell = params['ell']
    exp = params['exp']

    S = ProjectiveSpace(1, GF(D)).rational_points()

    P = choice(S)
    C2, D2 = P.dehomogenize(0), P.dehomogenize(1)
    mewtwo = StrongApproximationHeuristic(D, C2, D2, [(ell, exp)])
    while mewtwo is None:
        P = choice(S)
        C2, D2 = P.dehomogenize(0), P.dehomogenize(1)
        mewtwo = StrongApproximationHeuristic(D, C2, D2, [(ell, exp)])

    C0, D0 = EichlerModConstraint(mewtwo, I)

    D2_prime = Zmod(D)
    C2_prime = mod(-D2_prime * C2 * inverse_mod(D2, D), D)

    C1 = CRT([C0, C2_prime], [N, D])
    D1 = CRT([D0, D2_prime], [N, D])

    mewone = StrongApproximationHeuristic(N * D, C1, D1, [(ell, exp)])

    return mewone, mewtwo

# Finds an element of wanted norm in the Eichler suborder ZZ + D*I
def EichlerSuborderNormEquation(params, D, I):
    N = I.norm()

    assert gcd(N, D) == 1, 'N and D are not coprime'

    mewone, mewtwo = EichlerSuborderNormEquation_helper(params, D, I, N)
    while mewone is None:
        mewone, mewtwo = EichlerSuborderNormEquation_helper(params, D, I, N)

    return mewtwo * mewone

# Computes a set of 3 elements in the suborder ZZ + D*O, which generate it
def SmoothGen(params, O, D):
    B = params['B']
    O0 = params['O0']

    assert O.discriminant() == B.discriminant(), 'O is not a maximal order'
    assert D.is_prime(), 'D is not prime'

    L = set()
    I0 = ConnectingIdeal(params, O0, O)

    # Missing description of ZZ + D * O in target_order = ZZ + D * O, wait for Dr. Leroux's answer

    found_fam = False
    generating_fam = []
    while not found_fam:
        I, _, alpha = EquivalentPrimeIdealHeuristic(I0, random_elements=True)  # Same as RandomEquivalentPrimeIdeal
        while I is None:
            I, _, alpha = EquivalentPrimeIdealHeuristic(I0, random_elements=True)  # Same as RandomEquivalentPrimeIdeal

        theta = EichlerSuborderNormEquation(D, I)

        L.add(alpha * theta * alpha^(-1))

        if len(L) >= 3:
            for trial_generating_fam in combinations(L, 3):
                if target_order == B.ideal([1] + trial_generating_fam):
                    found_fam = True
                    generating_fam = trial_generating_fam
                    break

    return generating_fam


# Takes the ideal representation of an isogeny and converts it to suborder representation
def IdealToSuborder(params, I):
    B = params['B']
    O0 = params['O0']

    D = I.norm()
    O = I.left_order()
    O_prime = I.right_order()

    generating_fam = SmoothGen(B, O0, O, D)
    isogenies = []
    for theta_i in generating_fam:
        phi_i = IdealToIsogenyFromKLPT(O_prime * theta_i)
        isogenies.append(phi_i)

    return O, isogenies


def IdealSuborderNormEquation_helper(params, D, I, N, J, N_prime):
    ell = params['ell']
    exp = params['exp']

    S = ProjectiveSpace(1, GF(D)).rational_points()

    P = choice(S)
    C2, D2 = P.dehomogenize(0), P.dehomogenize(1)
    mewtwo = StrongApproximationHeuristic(D, C2, D2, [(ell, exp)])
    while mewtwo is None or gcd(mewtwo.norm(), N_prime) == 1:
        P = choice(S)
        C2, D2 = P.dehomogenize(0), P.dehomogenize(1)
        mewtwo = StrongApproximationHeuristic(D, C2, D2, [(ell, exp)])

    C0, D0 = EichlerModConstraint(mewtwo, I)
    C3, D3 = IdealModConstraint(mewtwo, J)

    D2_prime = Zmod(D).random_element()
    C2_prime = mod(-D2_prime * C2 * inverse_mod(D2, D), D)

    C1 = CRT([C0, C2_prime, C3], [N, D, N_prime])
    D1 = CRT([D0, D2_prime, D3], [N, D, N_prime])

    mewone = StrongApproximationHeuristic(N * D * N_prime, C1, D1, [(ell, exp)])

    return mewone, mewtwo

# Finds an element of wanted norm in some intersection of ideals (ZZ + D*I) `intersect` J
def IdealSuborderNormEquation(params, D, I, J):
    N = I.norm()
    N_prime = J.norm()

    assert gcd(N, N_prime) == 1 and gcd(N, D) == 1 and gcd(N_prime, D) == 1, "N, N' and D are not coprime"

    mewone, mewtwo = IdealSuborderNormEquation_helper(params, D, I, N, J, N_prime)
    while mewone is None:
        mewone, mewtwo = IdealSuborderNormEquation_helper(params, D, I, N, J, N_prime)

    return mewtwo * mewone


def powerset(iter, include_emptyset=True):
    lst = list(iter)
    length = len(lst)

    return chain.from_iterable(
        combinations(lst, size)
        for size in range(0 if include_emptyset else 1, length + 1)
    )

# Checks that the trace of the given endomorphisms all match their corresponding isogeny
def CheckTrace(params, M, E, isogenies, generating_fam):
    p = params['p']

    assert len(isogenies) == len(generating_fam), 'The number of isogenies and elements in the generating family should be the same'

    # Is E2 guaranteed to have 2 generators (i.e. it is the product of 2 cyclic groups)?
    P, Q = [((p + 1) / M) * G for G in E.gens()]  # From https://github.com/jack4818/Castryck-Decru-SageMath

    for I in powerset(range(len(isogenies)), include_emptyset=False):
        theta_I = generating_fam[0]
        phi_I = isogenies[0]
        for i in I[1:]:
            theta_I *= generating_fam[i]
            phi_I *= isogenies[i]

        for R in [P, Q]:
            phi_I_hat = phi_I.dual()
            if phi_I(R) + phi_I_hat(R) != theta_I.trace() * R:
                return False

    return True

# Given a suborder representation, runs public-key validation on it
def SuborderVerification(params, M, x, pi):
    p = params['p']
    O0 = params['O0']
    D, E1, E2 = x
    O, isogenies = pi

    if O.discriminant() != p:
        return False

    generating_fam = SmoothGen(params, O, D)
    J = ConnectingIdeal(params, O0, O)
    L = EquivalentPrimeIdealHeuristic(J)
    while L is None:
        L = EquivalentPrimeIdealHeuristic(J)

    psi = IdealToIsogenyFromKLPT(L)
    E1_prime = psi.codomain()
    if E1.j_invariant() != E1_prime.j_invariant() or E1.j_invariant() != E1_prime.j_invariant()^p:
        return False

    for phi_i, theta_i in zip(isogenies, generating_fam):
        assert phi_i.degree() == theta_i.norm(), 'The degree of the isogeny φᵢ should be the same as the norm of the element of the generating family θᵢ'

        F_i = phi_i.codomain()

        if F_i.j_invariant() != E2.j_invariant():
            return False

    return CheckTrace(params, M, E2, isogenies, generating_fam)


def find_large_enough_M(E, generating_fam):
    n = len(generating_fam)
    target_M = max([2 * sqrt(gen.norm()^n) for gen in generating_fam])

    m = 1
    M = E.count_points(m)[-1]
    while M <= target_M:
        m += 1
        M = E.count_points(m)[-1]

    return M

def multidimensional_discrete_log(generators, target):
    raise NotImplementedError()

# Given the other party's public parameters/suborder representation
# as well as your own, compute the kernel of the resulting isogeny composition
def SuborderEvaluation(params, E1, E2, pi, D, J):
    p = params['p']
    B = params['B']
    O0 = params['O0']
    O, isogenies = pi

    if J.left_order() != O:
        return None

    generating_fam = SmoothGen(B, O0, O, D)

    M = find_large_enough_M(E1, generating_fam)
    if not SuborderVerification(B, O0, M, (D, E1, E2), pi):
        return None

    L = ConnectingIdeal(B, O0, O)
    I, _, alpha = EquivalentPrimeIdealHeuristic(L, random_elements=True)  # Same as RandomEquivalentPrimeIdeal()?
    while I is None:
        I, _, alpha = EquivalentPrimeIdealHeuristic(L, random_elements=True)

    beta = IdealSuborderNormEquation(D, I, alpha^(-1) * J * alpha)

    generating_set = [prod(subset, B.one()) for subset in powerset(generating_fam)]
    coeffs = multidimensional_discrete_log(generating_set, alpha * beta * alpha.inverse())

    # Is E2 guaranteed to have 2 generators (i.e. it is the product of 2 cyclic groups)?
    P, Q = [((p + 1) / J.norm()) * G for G in E2.gens()]  # From https://github.com/jack4818/Castryck-Decru-SageMath

    R = S = E2(0)
    isogeny_compositions = [prod(subset, E2.isogeny(E2(0))) for subset in powerset(isogenies)]
    for coeff, isogeny in zip(coeffs, isogeny_compositions):
        R += coeff * isogeny(P)
        S += coeff * isogeny(Q)

    if S == 0:
        return Q

    a = DLP(R, S)

    return P - a * Q


# Generate a secret ideal representing an isogeny, and its public suborder representation
def KeyGeneration(params):
    B = params['B']
    O0 = params['O0']

    I = randomideal(O0)
    D = I.norm()
    while not D.is_prime():
        I = randomideal(O0)
        D = I.norm()

    pi = IdealToSuborder(params, I)
    _, isogenies = pi

    E = isogenies[0].domain()
    assert all([isogeny_i.domain() == E for isogeny_i in isogenies]), 'All isogenies should map from the curve E'

    return (E, pi), (I, D)


# Do the complete pSIDH key exchange
def KeyExchange(params, I, D_prime, E_prime, pi):
    p = params['p']
    B = params['B']
    E0 = params['E0']
    D = I.norm()

    assert D != D_prime, "D and D' should be different"

    O, isogenies = pi

    generating_fam = SmoothGen(params, O, D_prime)

    M = find_large_enough_M(E_prime, generating_fam)
    if not SuborderVerification(params, M, (D_prime, E0, E_prime), pi):
        return None

    J = O * 1

    theta = IdealSuborderNormEquation(D_prime, J, I)

    B = p^2 * D^9 * J.norm()^2
    T = B + 1
    min_smoothness_bound = factor(T)[-1][0]
    for potential_T in range(B + 2, 2 * B):
        smoothness_bound = factor(potential_T)[-1][0]
        if smoothness_bound < min_smoothness_bound:
            T = potential_T
            min_smoothness_bound = smoothness_bound
    T_facts = factor(T)

    G = set()
    for prime, multiplicity in T_facts:
        # Missing alpha in J_i = O * (alpha.inverse() * theta * alpha).conjugate() + O * prime^multiplicity, wait for Dr. Leroux's answer
        generating_point = SuborderEvaluation(params, E0, E_prime, pi, D_prime, J_i)
        if generating_point is None:
            return None
        G.add(generating_point)

    psi = E_prime.isogeny(list(G))

    return psi.codomain().j_invariant()
