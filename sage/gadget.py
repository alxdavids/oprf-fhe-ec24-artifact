"""
Gadget Matrices.
"""
from sage.all import ZZ, matrix, identity_matrix, vector, PolynomialRing


def gadget_matrix(n, B, ell, R=ZZ):
    Id = identity_matrix(R, n)
    g = matrix(R, 1, ell, [B**i for i in range(ell)])
    return Id.tensor_product(g)


def decompose_lwe(v, B, ell, R=ZZ):
    """

    EXAMPLE::

        >>> from sage.all import *
        >>> A = random_matrix(GF(127), 3, 4)
        >>> a = decompose_lwe(A, 2, 7)
        >>> G = gadget_matrix(4, 2, 7)
        >>> a*G.T == A
        True

    """
    try:
        _ = v[0, 0]  # is this a matrix?
        is_matrix = True
    except TypeError:
        is_matrix = False

    if is_matrix:
        return matrix(R, [decompose_lwe(v_, B, ell, R) for v_ in v.rows()])

    n = len(v)
    x = vector(ZZ, n * ell)
    for i in range(n):
        v_ = v[i].lift_centered()
        if v_ < 0:
            sgn = -1
            v_ = -v_
        else:
            sgn = 1
        for j in range(ell):
            x[i * ell + j] = sgn * (v_ % B)
            v_ = v_ // B
    return x


def decompose_rlwe(v, B, ell, d, R=PolynomialRing(ZZ, "x")):
    """

    EXAMPLE::

        >>> from sage.all import *
        >>> P = PolynomialRing(GF(127), "x")
        >>> A = matrix(P, 3, 4, [P.random_element(degree=3) for _ in range(3*4)])
        >>> a = decompose_rlwe(A, 2, 7, 4)
        >>> G = gadget_matrix(4, 2, 7)
        >>> a*G.T == A
        True


    """
    try:
        _ = v[0, 0]  # is this a matrix?
        is_matrix = True
    except TypeError:
        is_matrix = False

    if is_matrix:
        return matrix(R, [decompose_rlwe(v_, B, ell, d=d, R=R) for v_ in v.rows()])

    X = R.gen()
    n = len(v)
    w = vector(R, n * ell)
    for i in range(n):
        for j in range(d):
            v_ = v[i][j].lift_centered()
            if v_ < 0:
                sgn = -1
                v_ = -v_
            else:
                sgn = 1
            for k in range(ell):
                w[i * ell + k] += (sgn * (v_ % B)) * X**j
                v_ = v_ // B
    return w
