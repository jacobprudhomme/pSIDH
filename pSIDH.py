from itertools import chain, combinations
from sage.all import CRT, DLP, GF, ProjectiveSpace, ZZ, Zmod, choice, factor, gcd, inverse_mod, mod, prod

from deuring import randomideal
from sqisign.deuring import IdealToIsogenyFromKLPT
from sqisign.KLPT import EichlerModConstraint, EquivalentPrimeIdealHeuristic, IdealModConstraint, StrongApproximationHeuristic


def ConnectingIdeal(params, O1, O2):
    B = params['B']

    assert O1.discriminant() == B.discriminant(), 'O1 is not a maximal ideal'
    assert O2.discriminant() == B.discriminant(), 'O2 is not a maximal ideal'

    N = O1.intersection(O2).free_module().index_in(O1.free_module())  # Sage orders in quaternion algebras don't directly support index_in()
    I = N * O1 * O2

    assert I.left_order() == O1 and I.right_order() == O2, 'I is not an O1,O2-connecting ideal'

    return I


def EichlerSuborderNormEquation_helper(D, I, N):
    S = ProjectiveSpace(1, GF(D)).rational_points()

    P = choice(S)
    C2, D2 = P.dehomogenize(0), P.dehomogenize(1)
    mewtwo = StrongApproximationHeuristic(D, C2, D2)
    while mewtwo is None:
        P = choice(S)
        C2, D2 = P.dehomogenize(0), P.dehomogenize(1)
        mewtwo = StrongApproximationHeuristic(D, C2, D2)

    C0, D0 = EichlerModConstraint(mewtwo, I)

    D2_prime = Zmod(D)
    C2_prime = mod(-D2_prime * C2 * inverse_mod(D2, D), D)

    C1 = CRT([C0, C2_prime], [N, D])
    D1 = CRT([D0, D2_prime], [N, D])

    return StrongApproximationHeuristic(N * D, C1, D1), mewtwo

def EichlerSuborderNormEquation(D, I):
    N = I.norm()

    assert gcd(N, D) == 1, 'N and D are not coprime'

    mewone, mewtwo = EichlerSuborderNormEquation_helper(D, I, N)
    while mewone is None:
        mewone, mewtwo = EichlerSuborderNormEquation_helper(D, I, N)

    return mewtwo * mewone

def SmoothGen(params, O, D):
    B = params['B']
    O0 = params['O0']

    L = set()
    I0 = ConnectingIdeal(params, O0, O)

    target_order = ZZ + D * O

    found_fam = False
    generating_fam = []
    while not found_fam:
        I, _, alpha = EquivalentPrimeIdealHeuristic(I0, random_elements=True) # Same as RandomEquivalentPrimeIdeal()?
        while I is None:
            I, _, alpha = EquivalentPrimeIdealHeuristic(I0, random_elements=True) # Same as RandomEquivalentPrimeIdeal()?

        theta = EichlerSuborderNormEquation(D, I)

        L.add(alpha * theta * alpha^(-1))

        if len(L) >= 3:
            for trial_generating_fam in combinations(L, 3):
                if target_order == B.ideal([1] + trial_generating_fam):
                    found_fam = True
                    generating_fam = trial_generating_fam
                    break

    return generating_fam


def IdealToSuborder(params, I):
    B = params['B']
    O0 = params['O0']

    # assert I.quaternion_algebra() == QuaternionAlgebra(_, _, _)

    D = I.norm()
    O = I.left_order()
    O_prime = I.right_order()

    generating_fam = SmoothGen(B, O0, O, D)
    isogenies = []
    for theta_i in generating_fam:
        phi_i = IdealToIsogenyFromKLPT(O_prime * theta_i)
        isogenies.append(phi_i)

    return O, isogenies  # Wait, these are isogenies... Shouldn't they be suborders? And later compute isogenies using Velu's formulas?


def IdealSuborderNormEquation_helper(D, I, N, J, N_prime):
    S = ProjectiveSpace(1, GF(D)).rational_points()

    P = choice(S)
    C2, D2 = P.dehomogenize(0), P.dehomogenize(1)
    mewtwo = StrongApproximationHeuristic(D, C2, D2)
    while mewtwo is None or gcd(mewtwo.norm(), N_prime) == 1:
        P = choice(S)
        C2, D2 = P.dehomogenize(0), P.dehomogenize(1)
        mewtwo = StrongApproximationHeuristic(D, C2, D2)

    C0, D0 = EichlerModConstraint(mewtwo, I)
    C3, D3 = IdealModConstraint(mewtwo, J)

    D2_prime = Zmod(D).random_element()
    C2_prime = mod(-D2_prime * C2 * inverse_mod(D2, D), D)

    C1 = CRT([C0, C2_prime, C3], [N, D, N_prime])
    D1 = CRT([D0, D2_prime, D3], [N, D, N_prime])

    return StrongApproximationHeuristic(N * D * N_prime, C1, D1), mewtwo

def IdealSuborderNormEquation(D, I, J):
    N = I.norm()
    N_prime = J.norm()

    assert gcd(N, N_prime) == 1 and gcd(N, D) == 1 and gcd(N_prime, D) == 1, "N, N' and D are not coprime"

    mewone, mewtwo = IdealSuborderNormEquation_helper(D, I, N, J, N_prime)
    while mewone is None:
        mewone, mewtwo = IdealSuborderNormEquation_helper(D, I, N, J, N_prime)

    return mewtwo * mewone


def powerset(iter, include_emptyset=True):
    lst = list(iter)
    length = len(lst)

    return chain.from_iterable(
        combinations(lst, size)
        for size in range(0 if include_emptyset else 1, length + 1)
    )

def CheckTrace(M, E, isogenies, generating_fam):
    assert len(isogenies) == len(generating_fam), 'The number of isogenies and elements in the generating family should be the same'

    P, Q = E(0).division_points(M).basis()  # TURN INTO GROUP???

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

def SuborderVerification(params, M, x, pi):
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

    return CheckTrace(M, E2, isogenies, generating_fam)


def multidimensional_discrete_log(generators, target):
    raise NotImplementedError()

def SuborderEvaluation(params, E1, E2, pi, D, J):
    p = params['p']
    B = params['B']
    O0 = params['O0']
    O, isogenies = pi

    if J.left_order() != O:
        return None

    if not SuborderVerification(B, O0, E1.base_extend(GF(p^m)).order(), (D, E1, E2), pi):  # WHAT IS m HERE??? ANTONIN
        return None

    generating_fam = SmoothGen(B, O0, O, D)
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


def KeyExchange(params, I, D_prime, E_prime, pi):
    p = params['p']
    B = params['B']
    E0 = params['E0']
    D = I.norm()

    assert D != D_prime, "D and D' should be different"

    O, isogenies = pi

    generating_fam = SmoothGen(params, O, D_prime)

    if not SuborderVerification(params, E_prime.base_extend(GF(p^m)).order(), (D_prime, E0, E_prime), pi):  # WHAT IS m HERE??? ANTONIN
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
        # J_i = O0 * ___ + O0 * prime^multiplicity
        generating_point = SuborderEvaluation(p, B, O0, E0, E_prime, pi, D_prime, J_i)
        if generating_point is None:
            return None
        G.add(generating_point)

    psi = E_prime.isogeny(list(G))

    return psi.codomain().j_invariant()
