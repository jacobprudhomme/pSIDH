from itertools import combinations
from sage.all import ZZ

from sqisign.deuring import IdealToIsogenyFromKLPT
from sqisign.KLPT import EquivalentPrimeIdealHeuristic, EquivalentRandomEichlerIdeal, RepresentIntegerHeuristic
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
