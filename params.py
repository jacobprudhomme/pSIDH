from sage.all import EllipticCurve, GF, QuaternionAlgebra, var


def complete_params(params):
    assert params['p'] % 4 == 3, 'Prime must be p ≡ 3 mod 4 for a variety of reasons (e.g. supersingularity of the elliptic curve)'

    B = QuaternionAlgebra()
    i, j, k = B.gens()

    x = var('x')
    Fp2 = GF(params['p']^2, modulus=(x^2 + 1))  # Quadratic field of definition for the elliptic curve

    params['B'] = B  # Parent quaternion algebra
    params['E0'] = EllipticCurve(Fp2, [1, 0])  # Starting curve, supersingular and with known endomorphism ring
    params['O0'] = B.quaternion_order([1, i, (i + j) / 2, (1 + k) / 2])  # Maximal order corresponding to endomorphism ring of E0


# Small 54-bit prime p, nice for testing
p_toy = {
    'p': 9568331647090687,
    'ell': 2,  # Global small prime parameter used in places like StrongApproximation
    'exp': 213,  # Large enough exponent to find elements of wanted norm in Ideal/EichlerSuborderNormEquation
}
complete_params(p_toy)

# Prime p given in the SQISign paper, since the pSIDH paper does not provide any example parameters
p6983 = {
    'p': 73743043621499797449074820543863456997944695372324032511999999999999999999999,
    'ell': 2,  # Global small prime parameter used in places like StrongApproximation
    'exp': 1000,  # Large enough exponent to find elements of wanted norm in Ideal/EichlerSuborderNormEquation
}
complete_params(p6983)
