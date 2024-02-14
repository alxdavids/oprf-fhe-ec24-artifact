"""
Circuit Private Bootstrapping.

[Kluczniak22] Kluczniak, K. (2022). Circuit privacy for FHEW/TFHE-style fully homomorphic
  encryption in practice. Cryptology ePrint Archive, Report 2022/1459.
  https://eprint.iacr.org/2022/1459

"""

from sage.all import (
    IntegerModRing,
    PolynomialRing,
    ZZ,
    identity_matrix,
    ceil,
    floor,
    log,
    matrix,
    randint,
    round,
    vector,
)

from sage.stats.distributions.discrete_gaussian_integer import (
    DiscreteGaussianDistributionIntegerSampler,
)
from sage.stats.distributions.discrete_gaussian_lattice import (
    DiscreteGaussianDistributionLatticeSampler,
)

from gadget import gadget_matrix
from tfhe import GSW, GGSW, BlindRotation, ModSwitchBootstrap
from compression import PackedBlindRotation


class GadgetPreimageSampler:
    """
    Samples `x` such that `G ⋅ x = u mod q` for x in `[0,B-1]^ℓ`, u in ZZ
    """

    def __init__(self, B=2, q=2**5, sigma=2.0):
        self.B = B
        self.q = q
        self.ell = ceil(log(q, B))
        self.sigma = sigma
        self.G = gadget_matrix(1, B, self.ell, ZZ)
        self.basis = matrix(self.gadget_basis())

    def gadget_basis(self):
        """
        Gadget basis from MP12
        """
        basis_vectors = []
        for i in range(self.ell - 1):
            basis_vectors += [vector([0] * i + [self.B, -1] + [0] * (self.ell - i - 2))]

        if log(self.q, self.B) in ZZ:
            basis_vectors += [vector([0] * (self.ell - 1) + [self.B])]
        else:
            basis_vectors += [self.q.digits(self.B)]

        return basis_vectors

    def __call__(self, v):
        """
        sage: gadget = GadgetPreimageSampler(q=3*2**4)
        sage: P.<x> = PolynomialRing(IntegerModRing(gadget.q))
        sage: u = x^2 + x + 1
        sage: v = 34*x^3 + 23*x^2 + 5*x + 10
        sage: G = identity_matrix(ZZ,2).tensor_product(gadget.G)
        sage: G * gadget([u,v]).change_ring(P)
        (x^2 + x + 1, 34*x^3 + 23*x^2 + 5*x + 10)

        """
        try:
            _ = v[0, 0]  # is this a matrix?
            try:
                _ = v[0, 0].coefficients()
                R = PolynomialRing(ZZ, "x")
            except (AttributeError, IndexError):
                R = ZZ
            is_matrix = True
        except TypeError:
            try:
                _ = v[0].coefficients()
                R = PolynomialRing(ZZ, "X")
            except (AttributeError, IndexError):
                R = ZZ
            is_matrix = False

        if is_matrix:
            return matrix(R, [self.__call__(v_) for v_ in v.rows()])

        if R == ZZ:
            preimages = []
            for vi in v:
                coset = vector(ZZ, ZZ(vi).digits(self.B))
                if len(coset) != self.ell:
                    coset = vector(ZZ, list(coset) + [0] * (self.ell - len(coset)))
                D = DiscreteGaussianDistributionLatticeSampler(
                    self.basis, self.sigma, -1 * coset
                )
                preimages += list(D() + coset)
            return vector(ZZ, preimages)

        else:
            # the base ring of u is ZZ[X]
            X = R.gen()
            preimages = []
            for vi in v:
                coeffs = vector(ZZ, vi.coefficients(sparse=False))
                deg = len(coeffs)
                zz_preimage = list(self.__call__(coeffs))
                rearranged = matrix(
                    deg, self.ell, zz_preimage
                ).T  # a row is coeffs of a polynomial
                this_preimage = rearranged * vector([X**i for i in range(deg)])
                preimages += list(this_preimage)
            return vector(R, preimages)


class CPGSW(GSW):
    """
    Circuit-private GSW.
    """

    def __init__(
        self,
        n,
        q,
        chi_s="binary",
        chi_e=2.0,
        B=2,
        p=2,
        s=None,
        sigma_x=1.0,
        force_delta=False,
    ):
        super().__init__(
            n=n, q=q, chi_s=chi_s, chi_e=chi_e, B=B, p=p, s=s, force_delta=force_delta
        )
        self.gadget = GadgetPreimageSampler(B, q, sigma_x)
        self.sigma_x = sigma_x

    def mul(self, C0, C1):
        """
        Multiply two ciphertexts.

        :param C0: GSW ciphertext.
        :param C1: GSW ciphertext.

        EXAMPLE::

            sage: gsw = CPGSW(8, 2**12, B=4, p=4)
            sage: C0 = gsw(2)
            sage: C1 = gsw(3)
            sage: C = gsw.mul(C0, C1)
            sage: gsw.decrypt(C)
            2
            sage: gsw.decrypt(gsw.mul(C1,C0))
            2
            sage: C = gsw.mul(C1, C0+C1)
            sage: gsw.decrypt(C)
            3

        """
        C0 = matrix(C0.base_ring(), C0.nrows(), C0.ncols(), C0.list())
        C1 = self.gadget(C1)
        C = C1 * C0
        return C


class CPGGSW(GGSW):
    """
    Circuit-private GGSW.
    """

    def __init__(self, d, q, k, chi_s="binary", chi_e=2.0, B=2, p=2, s=None, sigma_x=1.0):
        super().__init__(d, q, k, chi_s, chi_e, B, p, s)
        self.gadget = GadgetPreimageSampler(B, q, sigma_x)
        self.sigma_x = sigma_x

    def mul(self, C0, C1):
        """
        EXAMPLE::

            sage: from tfhe import *
            sage: ggsw = GGSW(4, 2**10, 1, B=4, p=4)
            sage: C0 = ggsw([1,1,0,0])
            sage: C1 = ggsw([3,1,0,0])
            sage: ggsw.decrypt(C0)
            x + 1
            sage: ggsw.decrypt(C1)
            x + 3
            sage: C = ggsw.mul(C0, C1)
            sage: ggsw.decrypt(C)
            x^2 + 3

            sage: c1 = C1[-1] # GLWE ciphertext
            sage: GLWE.decrypt(ggsw, ggsw.mul(C0, c1))
            x^2 + 3

        """
        C0 = matrix(C0.base_ring(), C0.nrows(), C0.ncols(), C0.list())
        C = (self.gadget(C1).change_ring(C0.base_ring()) * C0) % self.phi
        return C


class CPBlindRotation(BlindRotation, CPGGSW):
    """
    Perform blind rotation with circuit privacy.

    EXAMPLE::

        sage: from tfhe import *
        sage: cpbr = CPBlindRotation(LWE(n=2, q=2^21, p=4), d=32, k=2)
        sage: c1 = cpbr(cpbr.lwe(1))
        sage: GLWE.decrypt(cpbr, c1)[0]
        1
        sage: c0 = cpbr(cpbr.lwe(0))
        sage: GLWE.decrypt(cpbr, c0)[0]
        0

        sage: cpbr = CPBlindRotation(LWE(n=10, q=2^20, p=4), d=32, k=2)
        sage: all([GLWE.decrypt(cpbr, cpbr(cpbr.lwe(1)))[0] == 1 for _ in range(3)])
        True

        sage: all([GLWE.decrypt(cpbr, cpbr(cpbr.lwe(0)))[0] == 0 for _ in range(3)])
        True

    """

    def __init__(
        self, lwe, d=1024, k=1, chi_s="binary", chi_e=2.0, B=2, p=2, s=None, sigma_x=1.0
    ):
        super().__init__(lwe, d, k, chi_s, chi_e, B, p, s)
        self.gadget = GadgetPreimageSampler(B, lwe.q, sigma_x)
        self.sigma_x = sigma_x


class CPModSwitchBootstrap(ModSwitchBootstrap, CPBlindRotation):
    def __init__(
        self,
        lwe,
        d=1024,
        k=1,
        chi_s="binary",
        chi_e=2.0,
        B=2,
        p_in=2,
        p_out=3,
        s=None,
        sigma_x=1.0,
        sigma_r=1.0,
        sigma_rand=1.0,
        LR=2,
    ):
        super().__init__(lwe, d, k, chi_s, chi_e, B, p_in, p_out, s)
        self.gadget = GadgetPreimageSampler(B, lwe.q, sigma_x)
        self.sigma_x = sigma_x

        self.lwe_ext = self.lwe.copy(n=k * d, p=p_out, chi_e=sigma_r, s=self.glwesec)
        self.LR = LR
        V = []
        lR = ceil(log(self.q, LR))
        for _ in range(lR):
            V.append(self.lwe_ext(0))
        self.V = matrix(self.lwe_ext.Rq, V)
        self.sigma_rand = sigma_rand
        self.lwe_o = self.lwe_ext

    def __call__(self, c):
        """
        Switch from mod 2 to mod 3 with circuit privacy (no keyswitching necessary)

        :param c: LWE ciphertext with plaintext modulus 2

        EXAMPLE::

            sage: from tfhe import *
            sage: cpmsbs = CPModSwitchBootstrap(LWE(n=20, q=2^10, p=2), d=64, k=1)
            sage: ct0 = cpmsbs(cpmsbs.lwe(0))
            sage: cpmsbs.lwe_o.decrypt(ct0)
            0

            sage: ct1 = cpmsbs(cpmsbs.lwe(1))
            sage: cpmsbs.lwe_o.decrypt(ct1)
            1

        """
        delta = self.q // self.p_out
        N = self.d
        v_modswitch = [delta] * (N // 2) + [-delta] * (N // 2)
        c = CPBlindRotation.__call__(self, c, v_modswitch)
        cext = self.sample_extract(c)

        gadget_rand = GadgetPreimageSampler(self.LR, self.q, self.sigma_rand)
        vecr = gadget_rand([0])
        Dr = DiscreteGaussianDistributionIntegerSampler(self.sigma_rand)
        Dy = DiscreteGaussianDistributionIntegerSampler(self.sigma_x)
        crand = vecr * self.V + vector([0] * (self.n * self.d) + [Dr()])
        cout = cext + crand + vector([0] * (self.n * self.d) + [Dy() - delta])

        return cout


class CPPackedBlindRotation(PackedBlindRotation, CPGGSW):
    """
    Perform blind rotation with circuit privacy with compressed
    keys. Only works for RLWE at the moment.

    EXAMPLE::

        sage: from tfhe import *
        sage: cpbr = CPPackedBlindRotation(LWE(n=2, q=2^20-1, p=4), d=32)
        sage: c1 = cpbr(cpbr.lwe(1))
        sage: GLWE.decrypt(cpbr, c1)[0]
        1
        sage: c0 = cpbr(cpbr.lwe(0))
        sage: GLWE.decrypt(cpbr, c0)[0]
        0

        sage: cpbr = CPBlindRotation(LWE(n=10, q=2^20, p=4), d=32)
        sage: all([GLWE.decrypt(cpbr, cpbr(cpbr.lwe(1)))[0] == 1 for _ in range(3)])
        True

        sage: all([GLWE.decrypt(cpbr, cpbr(cpbr.lwe(0)))[0] == 0 for _ in range(3)])
        True

    """

    def __init__(
        self, lwe, d=1024, k=1, chi_s="binary", chi_e=2.0, B=2, p=2, s=None, sigma_x=1.0
    ):
        super().__init__(lwe, d, k, chi_s, chi_e, B, p, s)
        self.gadget = GadgetPreimageSampler(B, lwe.q, sigma_x)
        self.sigma_x = sigma_x


class CPPackedModSwitchBootstrap(CPModSwitchBootstrap, CPPackedBlindRotation):
    pass
