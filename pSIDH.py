from itertools import combinations
from sage.all import CRT, ZZ, Zmod, gcd, inverse_mod, mod

from sqisign.deuring import IdealToIsogenyFromKLPT
from sqisign.KLPT import EichlerModConstraint, EquivalentPrimeIdealHeuristic, EquivalentRandomEichlerIdeal, IdealModConstraint, RepresentIntegerHeuristic, StrongApproximationHeuristic
from sqisign.setup import B, Bτ, O0, eτ, p, ω, l
from sqisign.utilities import inert_prime


def ConnectingIdeal(O1, O2):
    # Compute a random prime ≤ Bτ which is inert
    # in R[ω].
    # Note: this is the same as picking p ≡ 3 mod 4
    Nl = l^eτ

    # Stop infinite loops
    for _ in range(1000):
        Nτ = inert_prime(Bτ, -ZZ(ω^2))
        # We need the product to be large enough for
        # RepresentIntegerHeuristic.
        if Nτ * Nl > 2 * p:
            break

    # Compute an endomorphism γ of norm Nτ l^eτ
    # Nτ < Bτ
    γ = None

    # Stop infinite loops
    for _ in range(1000):
        γ = RepresentIntegerHeuristic(Nτ * Nl, parity=True)
        if γ is not None:
            break

    if γ is None:
        exit("Never found an alg element with norm (Nτ * Nl), Exiting...")

    I = O1 * γ + O2 * Nτ
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
        I = EquivalentPrimeIdealHeuristic(I0) # Same as RandomEquivalentPrimeIdeal()
        while I is None:
            I = EquivalentPrimeIdealHeuristic(I0) # Same as RandomEquivalentPrimeIdeal()

        alpha = I / I0
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
