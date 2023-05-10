from sage.all import ZZ

from sqisign.KLPT import RepresentIntegerHeuristic
from sqisign.setup import Bτ, eτ, p, ω, l
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
