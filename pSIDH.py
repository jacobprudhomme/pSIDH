from itertools import chain, combinations
from sage.all import CRT, DLP, GF, ZZ, Zmod, factor, gcd, inverse_mod, mod

from deuring import randomideal
from sqisign.deuring import IdealToIsogenyFromKLPT
from sqisign.KLPT import EichlerModConstraint, EquivalentPrimeIdealHeuristic, EquivalentRandomEichlerIdeal, IdealModConstraint, StrongApproximationHeuristic
from sqisign.setup import B, E0, O0, p


def ConnectingIdeal(B, O1, O2):
    assert O1.discriminant() == B.discriminant(), 'O1 is not a maximal ideal'
    assert O2.discriminant() == B.discriminant(), 'O2 is not a maximal ideal'

    N = O1.intersection(O2).free_module().index_in(O1.free_module())  # Sage orders in quaternion algebras don't directly support index_in()
    I = N * O1 * O2

    assert I.left_order() == O1 and I.right_order() == O2, 'I is not an O1,O2-connecting ideal'

    return I


def SmoothGen(O, D):
    # N = a subset of the natural numbers following certain constraints
    # O = maximal order
    # D = prime
    # Returns a generating family theta_1, theta_2, theta_3 for Z + DO, where n(theta_j) is in N

    L = set()
    I0 = ConnectingIdeal(O0, O)

    target_order = ZZ + D * O

    found_fam = False
    generating_fam = []
    while not found_fam:
        I, _, alpha = EquivalentPrimeIdealHeuristic(I0, random_elements=True) # Same as RandomEquivalentPrimeIdeal()?
        while I is None:
            I, _, alpha = EquivalentPrimeIdealHeuristic(I0, random_elements=True) # Same as RandomEquivalentPrimeIdeal()?

        theta = EquivalentRandomEichlerIdeal(I, D) # Same as EichlerSuborderNormEquation()

        L.add(alpha * theta * alpha^(-1))

        if len(L) >= 3:
            for trial_generating_fam in combinations(L, 3):
                if target_order == B.ideal([1] + trial_generating_fam):
                    found_fam = True
                    generating_fam = trial_generating_fam
                    break

    return generating_fam


def IdealToSuborder(I):
    # assert(I.quaternion_algebra() == QuaternionAlgebra(_, _, _))

    D = I.norm()
    O = I.left_order()
    O_prime = I.right_order()

    generating_fam = SmoothGen(O, D)
    isogenies = []
    for theta_i in generating_fam:
        phi_i = IdealToIsogenyFromKLPT(O_prime * theta_i)
        isogenies.append(phi_i)

    return O, isogenies  # Wait, these are isogenies... Shouldn't they be suborders? And later compute isogenies using Velu's formulas?


def IdealSuborderNormEquation_helper(D, I, N, J, N_prime):
    # SELECT A RANDOM CLASS???
    mewtwo = StrongApproximationHeuristic(D, C2, D2)
    while mewtwo is None or gcd(mewtwo.norm(), N_prime) == 1:
        # SELECT A RANDOM CLASS???
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

    assert(gcd(N, N_prime) == 1 and gcd(N, D) == 1 and gcd(N_prime, D) == 1)

    mewone, mewtwo = IdealSuborderNormEquation_helper(D, I, N, J, N_prime)
    while mewone is None:
        mewone, mewtwo = IdealSuborderNormEquation_helper(D, I, N, J, N_prime)

    return mewone * mewtwo


def all_proper_nonempty_subsets(iter):
    lst = list(iter)
    length = len(lst)

    return chain.from_iterable(
        combinations(lst, size)
        for size in range(1, length)
    )

def CheckTrace(M, E, isogenies, generating_fam):
    assert(len(isogenies) == len(generating_fam))

    P, Q = E(0).division_points(M).basis()  # TURN INTO GROUP???

    for I in all_proper_nonempty_subsets(range(len(isogenies))):  # IS IT REALLY JUST PROPER SUBSETS??? ANTONIN
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

def SuborderVerification(M, x, pi):
    D, E1, E2 = x
    O, isogenies = pi

    if O.discriminant() != p:
        return False

    generating_fam = SmoothGen(O, D)
    J = ConnectingIdeal(O0, O)
    L = EquivalentPrimeIdealHeuristic(J)
    while L is None:
        L = EquivalentPrimeIdealHeuristic(J)

    psi = IdealToIsogenyFromKLPT(L)
    E1_prime = psi.codomain()
    if E1.j_invariant() != E1_prime.j_invariant() or E1.j_invariant() != E1_prime.j_invariant()^p:
        return False

    for phi_i, theta_i in zip(isogenies, generating_fam):
        assert(phi_i.degree() == theta_i.norm())

        F_i = phi_i.codomain()

        if F_i.j_invariant() != E2.j_invariant():
            return False

    return CheckTrace(M, E2, isogenies, generating_fam)

def SuborderEvaluation(E1, E2, pi, D, J):
    O, isogenies = pi

    if J.left_order() != O:
        return None

    if not SuborderVerification(E1.base_extend(GF(p^m)).order(), (D, E1, E2), pi):  # WHAT IS m HERE??? ANTONIN
        return None

    generating_fam = SmoothGen(O, D)
    L = ConnectingIdeal(O0, O)
    I, _, alpha = EquivalentPrimeIdealHeuristic(L, random_elements=True)  # Same as RandomEquivalentPrimeIdeal()?
    while I is None:
        I, _, alpha = EquivalentPrimeIdealHeuristic(L, random_elements=True)

    beta = IdealSuborderNormEquation(D, I, alpha^(-1) * J * alpha)

    # EXPRESS SOMETHING AS LINEAR COMBINATION OF GENERATING FAMILY??? ANTONIN

    P, Q = E2(0).division_points(J.norm()).basis()  # TURN INTO GROUP???

    # EXPRESS SOMETHING AS LINEAR COMBINATION OF ISOGENIES OVER BASIS??? ANTONIN

    if S == 0:
        return Q

    a = DLP(R, S)

    return P - a * Q


def KeyGeneration():
    I = randomideal(O0)
    D = I.norm()
    while not D.is_prime():
        I = randomideal(O0)
        D = I.norm()

    pi = IdealToSuborder(I)
    _, isogenies = pi

    E = isogenies[0].domain()
    assert all([isogeny_i.domain() == E for isogeny_i in isogenies])

    return (E, pi), (I, D)


def KeyExchange(I, D_prime, E_prime, pi):
    D = I.norm()

    assert(D != D_prime)

    O0, isogenies = pi

    generating_fam = SmoothGen(O0, D_prime)

    if not SuborderVerification(E_prime.base_extend(GF(p^m)).order(), (D_prime, E0, E_prime), pi):  # WHAT IS m HERE??? ANTONIN
        return None

    J = O0 * 1

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
        generating_point = SuborderEvaluation(E0, E_prime, pi, D_prime, J_i)  # THIS CAN BE None, HOW SHOULD THIS BE HANDLED??? ANTONIN
        if generating_point is None:
            return None
        G.add(generating_point)

    psi = E_prime.isogeny(list(G))

    return psi.codomain().j_invariant()
