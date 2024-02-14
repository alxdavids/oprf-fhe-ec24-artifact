"""
Toy Implementation of TFHE.

LITERATURE:

[CGGI20] Chillotti, I., Gama, N., Georgieva, M., & Izabachène, M. (2020). TFHE: fast fully
  homomorphic encryption over the torus. Journal of Cryptology, 33(1), 34–91.
  http://dx.doi.org/10.1007/s00145-019-09319-x

[Joye21] Joye, M. (2021). Guide to fully homomorphic encryption over the [discretized]
  torus. Cryptology ePrint Archive, Report 2021/1402. https://eprint.iacr.org/2021/1402

"""

from sage.all import (
    IntegerModRing,
    PolynomialRing,
    ZZ,
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

from gadget import gadget_matrix, decompose_lwe, decompose_rlwe


def apply_plaintext_matrix(A, other, balance=True, randomize=False):
    """
    Apply this matrix to vector of elements, considering this matrix over ZZ.

    :param other: some iterable.
    :param balance: consider elements in `(-q//2, …, q//2]`.
    :param randomize: randomize with ± when p=2

    EXAMPLE::

        sage: p = 2
        sage: A = random_matrix(GF(p), 5, 8)
        sage: lwe = LWE(10, 127, "binary", 3.0, p=p)
        sage: x = random_vector(GF(p), 8)
        sage: y = A*x
        sage: c0 = [lwe(x_) for x_ in x]
        sage: c1 = apply_plaintext_matrix(A, c0)
        sage: c2 = apply_plaintext_matrix(A, c0)
        sage: c3 = apply_plaintext_matrix(A, c0, randomize=True)
        sage: z = vector(GF(p), [lwe.decrypt(c_) for c_ in c1])
        sage: y == z
        True
        sage: z = vector(GF(p), [lwe.decrypt(c_) for c_ in c2])
        sage: y == z
        True
        sage: z = vector(GF(p), [lwe.decrypt(c_) for c_ in c3])
        sage: y == z
        True
        sage: c1 == c2
        True
        sage: c2 == c3
        False

    """

    res = [0] * A.nrows()

    if randomize is True and A.base_ring().characteristic() != 2:
        raise ValueError(f"Cannot randomize signs in {A.base_ring()}.")

    for i in range(A.nrows()):
        for j in range(A.ncols()):
            if balance:
                c = A[i, j].lift_centered()
            else:
                c = A[i, j].lift()
            if randomize:
                c = (-1) ** ZZ.random_element(2) * c
            res[i] += c * other[j]
    return tuple(res)


class LWE:
    """
    LWE Distribution.

    EXAMPLE::

        sage: lwe = LWE(10, 127, "binary", 2.0)
        sage: _ = lwe()

    """

    def lift_centered(self, e):
        """
        Lift to the Integers

        :param e: some element

        EXAMPLE::

            sage: lwe = LWE(10, 127, "binary", 2.0)
            sage: c = lwe()
            sage: cz = lwe.lift_centered(c)
            sage: max(cz) <= 127//2
            True
            sage: min(cz) < 0
            True
            sage: lwe.lift_centered(GF(127)(126))
            -1
        """
        try:
            return e.lift_centered()
        except AttributeError:
            return self.R([e_.lift_centered() for e_ in list(e)])

    def copy(self, **kwds):
        """
        EXAMPLE::

            sage: lwe = LWE(10, 127)
            sage: lwe.copy(n=11, q=128)
            LWE(n=11, q=128, χ_e=2, p=2)
        """
        return LWE(
            n=kwds.get("n", self.n),
            q=kwds.get("q", self.q),
            chi_s=kwds.get("chi_s", self.chi_s),
            chi_e=kwds.get("chi_e", self.chi_e),
            p=kwds.get("p", self.p),
            s=kwds.get("s", self._s[: self.n]),
        )

    def __repr__(self):
        return f"LWE(n={self.n}, q={self.q}, χ_e={self.e.sigma:.1}, p={self.p})"

    @staticmethod
    def normalize_distribution(D):
        """
        Turn user-friendly descriptions of distributions into objects we can use.

        :param D: "binary", "ternary" or a standard deviation

        EXAMPLE::

            sage: D = LWE.normalize_distribution("binary")
            sage: max([D() for _ in range(1000)]) == 1
            True
            sage: min([D() for _ in range(1000)]) == 0
            True

            sage: D = LWE.normalize_distribution("ternary")
            sage: max([D() for _ in range(1000)]) == 1
            True
            sage: min([D() for _ in range(1000)]) == -1
            True

            sage: D = LWE.normalize_distribution(3.0)
            sage: D
            Discrete Gaussian sampler over the Integers with sigma = 3.000000 and c = 0.0...


        """
        if D == "binary":
            return lambda: randint(0, 1)
        if D == "ternary":
            return lambda: randint(-1, 1)
        try:
            return DiscreteGaussianDistributionIntegerSampler(sigma=D)
        except TypeError:
            return D

    def __init__(self, n, q, chi_s="binary", chi_e=2.0, p=2, s=None):
        """

        :param n: LWE secret dimension
        :param q: Modulus ∈ ZZ and ≥ 2
        :param chi_s: Secret distribution
        :param chi_e: Error distribution
        :param p: Plaintext modulus
        :param s: Explicitly set a secret.

        EXAMPLE::

            sage: lwe = LWE(10, 127)
            sage: lwe = LWE(10, 2^7)
            sage: lwe = LWE(10, 17, s=[0])

        """

        if p > q:
            raise ValueError(
                f"Plaintext space {p} is too big relative to ciphertext modulus {q}."
            )

        self.n = n
        self.q = q
        self.Rq = IntegerModRing(q)
        self.R = ZZ
        chi_s = self.normalize_distribution(chi_s)
        self.chi_s = chi_s
        self.phi = q
        self.e = self.normalize_distribution(chi_e)
        self.chi_e = chi_e
        if s is None:
            self._s = vector([chi_s() for _ in range(n)] + [1])
        else:
            if s == 0:
                s = [0] * self.n
            self._s = vector(list(s) + [1])
        self.p = p
        self.delta = floor(q / float(p))

    def _random_scalar(self):
        """
        Return a random scalar.

        EXAMPLE::

            sage: lwe = LWE(10, 127)
            sage: lwe._random_scalar() in GF(127)
            True

        """
        return self.Rq.random_element()

    def a(self):
        """
        Sample `a` component of an LWE ciphertext.

        EXAMPLE::

            sage: lwe = LWE(10, 127)
            sage: lwe.a() in VectorSpace(GF(127), 10)
            True

        """
        return vector([self._random_scalar() for _ in range(self.n)])

    def __call__(self, m=0, raw=False):
        """
        Encrypt `m`.

        :param m: m ∈ ZZ_p
        :param raw: do not encode the value (i.e. m ∈ ZZ_q)

        """
        a = self.a()
        e = self.e()
        b = (a * self._s[: self.n]) % self.phi + e
        if m:
            if raw:
                b += m
            else:
                m = self.R(self.R(m) % self.p)
                b += self.delta * m

        ab = vector(self.Rq, self.n + 1, list(-a) + [b])
        return ab

    def decrypt(self, c, raw=False):
        """
        Decrypt ciphertext `c`.

        :param c: LWE ciphertext
        :param raw: Do not decode.

        EXAMPLE::

        sage: lwe = LWE(10, 127, "binary", 2.0)
        sage: c = lwe(m=1)
        sage: lwe.decrypt(c)
        1

        sage: lwe = LWE(10, 127, "binary", 2.0, p=7)
        sage: c = lwe(m=6)
        sage: lwe.decrypt(c)
        6
        sage: round(lwe.decrypt(c,raw=True) / lwe.delta) % lwe.p
        6

        """
        z = self.lift_centered(c * self._s)
        if raw:
            return z
        else:
            return round(z / self.delta) % self.p


class GLWE(LWE):
    """
    Generalised LWE Distribution.
    """

    def __init__(self, d, q, k=1, chi_s="binary", chi_e=2.0, p=2, s=None):
        """
        :param d: Ring dimension, not enforced to a power of two, but it should be
        :param q: Modulus ∈ ZZ and ≥ 2
        :param chi_s: Secret distribution
        :param chi_e: Error distribution
        :param p: Plaintext modulus ∈ ZZ
        :param s: Explicitly set a secret.

        EXAMPLE::

            sage: glwe = GLWE(8, 127, 3)
            sage: glwe = GLWE(8, 2^7, 3)
            sage: glwe = GLWE(8, 17, s=[0])

        """
        self.n = k
        self.d = d
        self.q = q
        self.Rq = PolynomialRing(IntegerModRing(q), "x")
        self.R = PolynomialRing(ZZ, "x")
        self.phi = self.R.gen() ** d + 1
        chi_s = self.normalize_distribution(chi_s)
        self.e = self.normalize_distribution(chi_e)

        if s is None:
            self._s = []
            for _ in range(self.n):
                self._s += [self.R([chi_s() for _ in range(self.d)])]
            self._s += [self.R(1)]
            self._s = vector(self._s)
        else:
            self._s = vector([self.R(s_) for s_ in s] + [self.R(1)])

        try:
            self.delta = floor(q / float(p))
            self.p = p
        except:
            self.p = 2
            self.delta = floor(q / 2)

    def __repr__(self):
        """
        EXAMPLE::

            sage: GLWE(16, 127, 2)
            GLWE(d=16, q=127, χ_e=2, p=2)

        """
        return f"GLWE(d={self.d}, q={self.q}, χ_e={self.e.sigma:.1}, p={self.p})"

    def _random_scalar(self):
        """
        Return a random scalar.

        EXAMPLE::

            sage: glwe = GLWE(8, 127, 2)
            sage: glwe._random_scalar() in PolynomialRing(GF(127), "x")
            True

        """
        return self.Rq.random_element(degree=self.d - 1)

    def lift_centered(self, e):
        return self.R([e_.lift_centered() for e_ in list(e)])

    def decrypt(self, c, raw=False):
        """
        Decrypt ciphertext `c`.

        :param c: Ciphertext
        :param raw: Do not decode.

        EXAMPLE::

            sage: glwe = GLWE(4, 127, 3)
            sage: c = glwe(m=[1,0,0,1])
            sage: glwe.decrypt(c)
            x^3 + 1
            sage: glwe.decrypt(c, raw=True) # random
            63*x^3 + 63

        """
        z = self.lift_centered(c * self._s % self.phi)
        if raw:
            return z
        z = z / self.delta
        return self.R([round(z_) % self.p for z_ in list(z)])


class KeySwitching:
    """
    Switch between LWE keys.
    """

    def __init__(self, lwe_out, lwe_in, B=2):
        """
        :param lwe_out: Output LWE instance
        :param lwe_in: Source LWE instance
        :param B: Decomposition base

        """
        if lwe_in.q != lwe_out.q:
            raise ValueError(f"Modulus mismatch: {lwe_in.q} != {lwe_out.q}")
        self.B = B
        self.ell = ceil(log(lwe_in.q, B))
        self.lwe_i = lwe_in
        self.lwe_o = lwe_out
        self.ksk = self.key_gen()

    def key_gen(self):
        ksk_list = [[[] for j in range(self.ell)] for i in range(self.lwe_i.n)]
        for i in range(self.lwe_i.n):
            for j in range(self.ell):
                ksk_list[i][j] = self.lwe_o(self.lwe_i._s[i] * self.B**j, raw=True)
            ksk_list[i] = tuple(ksk_list[i])
        return tuple(ksk_list)

    def __call__(self, c):
        """
        Switch ciphertext `c` to output LWE instance.

        EXAMPLE::

            sage: lwe0 = LWE(10, 2047, p=7)
            sage: lwe1 = LWE(9, 2047, p=7)
            sage: ks = KeySwitching(lwe_out=lwe1, lwe_in=lwe0)
            sage: c0 = lwe0(5)
            sage: lwe1.decrypt(ks(c0))
            5
            sage: c1 = lwe0(2)
            sage: lwe1.decrypt(ks(c1))
            2

        """
        a, b = c[:-1], c[-1]
        ctxt_out = vector([0] * self.lwe_o.n + [b])

        for i in range(self.lwe_i.n):
            try:
                ai = decompose_lwe([a[i]], self.B, self.ell)
            except:
                ai = decompose_rlwe([a[i]], self.B, self.ell, self.lwe_i.d)
            for j in range(self.ell):
                ctxt_out += ai[j] * self.ksk[i][j]

        return ctxt_out


class GSW(LWE):
    """
    GSW Distribution.

    """

    def __init__(
        self, n, q, chi_s="binary", chi_e=2.0, B=2, p=2, s=None, force_delta=False
    ):
        """
        :param n: LWE secret dimension
        :param q: Modulus ∈ ZZ and ≥ 2
        :param chi_s: Secret distribution
        :param chi_e: Error distribution
        :param p: Plaintext modulus
        :param s: Explicitly set a secret.
        :param force_delta: throw an error if G[-1, -1] ≠ δ

        EXAMPLE::

            sage: gsw = GSW(10, 128)

        """

        super().__init__(n=n, q=q, chi_s=chi_s, chi_e=chi_e, p=p, s=s)
        self.B = B
        self.ell = ceil(log(q, B))
        self.G = gadget_matrix(self.n + 1, self.B, self.ell, self.R)
        if force_delta and self.G[-1, -1] != self.delta:
            raise ValueError(f"G[-1,-1] = {self.G[-1,1]} ≠ {self.delta} = δ.")

    def base_scheme(self):
        """
        Return base LWE instance with matching parameters and secret.

        """
        lwe = LWE(self.n, self.q, p=self.p, s=self._s[:-1])
        lwe.e = self.e
        return lwe

    def __call__(self, m=0):
        """
        Encrypt `m`

        :param m: m ∈ ZZ_p

        EXAMPLE::

            sage: gsw = GSW(9, 256, B=4, p=4)
            sage: C0 = gsw(2)
            sage: C1 = gsw(3)

        """
        Z = []
        for _ in range((self.n + 1) * self.ell):
            Z.append(super().__call__())
        Z = matrix(Z)
        if m:
            m = self.R(self.R(m) % self.p)
            C = Z + m * self.G.T
        else:
            C = Z
        return C

    def decrypt(self, c):
        """
        Decrypt `c`.

        :param c: Ciphertext.

        EXAMPLE::

            sage: gsw = GSW(9, 256, B=4, p=4)
            sage: gsw.decrypt(gsw(2))
            2
            sage: gsw.decrypt(gsw(3))
            3

        """
        return super().decrypt(c[-1])

    def mul(self, C0, C1):
        """
        Multiply two ciphertext.

        :param C0: GSW ciphertext.
        :param C1: GSW ciphertext.

        EXAMPLE::

            sage: gsw = GSW(8, 1024, B=4, p=4)
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

            sage: c1 = C1[-1] # LWE ciphertext
            sage: LWE.decrypt(gsw, gsw.mul(C0, c1))
            2

        """
        C0 = matrix(C0.base_ring(), C0.nrows(), C0.ncols(), C0.list())
        C1 = decompose_lwe(C1, self.B, self.ell, self.R)
        C = C1 * C0
        return C

    def cmux(self, Cb, c0, c1):
        """
        Select `c0` or `c1` depending on bit encrypted under `Cb`.

        :param Cb: GSW ciphertext of selector bit.
        :param c0: LWE ciphertext.
        :param c1: LWE ciphertext.

        EXAMPLE::

            sage: gsw = GSW(4, 1024, chi_s="binary", B=4, p=4)
            sage: C0 = gsw(0)
            sage: C1 = gsw(1)
            sage: lwe = gsw.base_scheme()
            sage: c0 = lwe(0)
            sage: c1 = lwe(1)
            sage: lwe.decrypt(gsw.cmux(C0, c0, c1))
            0
            sage: lwe.decrypt(gsw.cmux(C1, c0, c1))
            1

        """

        return self.mul(Cb, c1 - c0) + c0


class GGSW(GLWE):
    def __init__(
        self, d, q, k, chi_s="binary", chi_e=2.0, B=2, p=2, s=None, force_delta=False
    ):
        """
        :param d: MLWE ring dimension
        :param q: Modulus ∈ ZZ and ≥ 2
        :param k: MLWE module rank
        :param chi_s: Secret distribution
        :param chi_e: Error distribution
        :param B: Decomposition base
        :param p: Plaintext modulus
        :param s: Explicitly set a secret.
        :param force_delta: throw an error if G[-1, -1] ≠ δ

        EXAMPLE::

            sage: gsw = GGSW(8, 128, 2)

        """
        super().__init__(d=d, q=q, k=k, chi_s=chi_s, chi_e=chi_e, p=p, s=s)
        self.B = B
        self.ell = ceil(log(q, B))
        self.G = gadget_matrix(self.n + 1, self.B, self.ell, self.R)

        if force_delta and self.G[-1, -1] != self.delta:
            raise ValueError(f"G[-1,-1] = {self.G[-1,1]} ≠ {self.delta} = δ.")

    def base_scheme(self):
        glwe = GLWE(self.d, self.q, k=self.n, p=self.p, s=self._s[: self.n])
        glwe.e = self.e
        return glwe

    def __call__(self, m=None):
        Z = []
        for _ in range((self.n + 1) * self.ell):
            Z.append(super().__call__())
        Z = matrix(self.Rq, Z)
        if m:
            m = self.R(self.R(m) % self.p)
            C = Z + (m * self.G.T) % self.phi
        else:
            C = Z
        return C

    def decrypt(self, c):
        return super().decrypt(c[-1])

    def mul(self, C0, C1):
        """
        EXAMPLE::

            sage: ggsw = GGSW(4, 2**10, 3, B=4, p=4)
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
        C1 = decompose_rlwe(C1, self.B, self.ell, self.d, self.R)
        C = (C1 * C0) % self.phi
        return C

    def cmux(self, Cb, c0, c1):
        """

        :param Cb:
        :param c0:
        :param c1:

        EXAMPLE::

            sage: ggsw = GGSW(6, 1024, k=2, B=4, p=4)
            sage: C0 = ggsw(0)
            sage: C1 = ggsw(1)
            sage: glwe = ggsw.base_scheme()
            sage: c2 = glwe(2)
            sage: c3 = glwe(3)
            sage: glwe.decrypt(ggsw.cmux(C0, c2, c3))
            2
            sage: glwe.decrypt(ggsw.cmux(C1,c2, c3))
            3
        """
        return self.mul(Cb, c1 - c0) + c0


class BlindRotation(GGSW):
    def __init__(self, lwe, d=1024, k=1, chi_s="binary", chi_e=2.0, B=2, p=2, s=None):
        """
        :param lwe: LWE instance
        :param d: MLWE ring dimension
        :param q: Modulus ∈ ZZ and ≥ 2
        :param k: MLWE module rank
        :param chi_s: Secret distribution
        :param chi_e: Error distribution
        :param B: Decomposition base
        :param p: Plaintext modulus
        :param s: Explicitly set a secret.

        """
        super().__init__(d, lwe.q, k, chi_s, chi_e, B, p, s)
        self.lwe = lwe
        self.brk = self.key_gen()

    def key_gen(self):
        brk_list = []
        for j in range(self.lwe.n):
            brk_list.append(super().__call__(self.lwe._s[j]))  # scalar -> constant coeff
        return brk_list

    def modulus_switch(self, ctxt):
        """

        EXAMPLE::

            sage: br = BlindRotation(LWE(n=10, q=2^10, p=2))
            sage: N = 2*br.d
            sage: ctxt = br.modulus_switch(br.lwe(1))
            sage: round((ctxt * br.lwe._s).lift() / (N/br.lwe.p))
            1
            sage: ctxt = br.modulus_switch(br.lwe(0))
            sage: round((ctxt * br.lwe._s).lift() / (N/br.lwe.p)) % br.lwe.p
            0

        """
        N = ZZ(self.d)
        q = ZZ(self.q)
        ctxt = ctxt.lift()
        ctxt = vector(ZZ, [round(2 * N / q * c_) for c_ in ctxt])
        ctxt = ctxt.change_ring(IntegerModRing(2 * N))
        return ctxt

    def test_polynomial(self):
        """
        Return the standard test polynomial

        :param fun:

        EXAMPLE::

            sage: br = BlindRotation(LWE(n=10, q=2^10, p=2), d=8)
            sage: br.test_polynomial()
            512*x^5 + 512*x^4 + 512*x^3 + 512*x^2

        """
        p = ZZ(self.p)
        N = ZZ(self.d)

        v = []
        for j in range(self.d):
            j = (j + self.d // 4) % self.d
            v.append(self.delta * (round((p * j) / (2 * N)) % p))
        return self.Rq(v) % self.phi

    def __call__(self, c, v=None, in_clear=False):
        """
        Perform blind rotation.

        :param ctxt:
        :param v: Custom test polynomial

        EXAMPLE::

            sage: br = BlindRotation(LWE(n=2, q=2^30, p=4), d=1024, k=2)
            sage: c1 = br(br.lwe(1))
            sage: GLWE.decrypt(br, c1)[0]
            1
            sage: c0 = br(br.lwe(0))
            sage: GLWE.decrypt(br, c0)[0]
            0

            sage: br = BlindRotation(LWE(n=10, q=2^20, p=4), d=32, k=2)
            sage: all([GLWE.decrypt(br, br(br.lwe(1)))[0] == 1 for _ in range(32)])
            True

            sage: all([GLWE.decrypt(br, br(br.lwe(0)))[0] == 0 for _ in range(32)])
            True


        """

        c = self.modulus_switch(c)
        a, b = c[: self.lwe.n], c[self.lwe.n]
        if v is None:
            v = self.test_polynomial()

        X = self.R.gen()
        acc = X ** (-b) * vector(self.Rq, [0] * self.n + [v]) % self.phi
        for j in range(self.lwe.n):
            if in_clear:
                if self.lwe._s[j] == 1:
                    acc = (X ** (-a[j]) * acc) % self.phi
            else:
                c0 = acc
                c1 = (X ** (-a[j]) * acc) % self.phi
                acc = self.cmux(self.brk[j], c0, c1)
        return acc


class ModSwitchBootstrap(BlindRotation):
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
        ks=False,
    ):
        super().__init__(lwe, d, k, chi_s, chi_e, B, p_in, s=None)
        self.p_out = p_out
        glwesec = []
        for i in range(k):
            rlwesec = self._s[i].coefficients(sparse=False)
            glwesec += rlwesec + [0] * (d - len(rlwesec))
        self.glwesec = glwesec
        self.ksbool = ks
        if ks:
            self.lwe_o = self.lwe.copy(p=p_out)
            self.ks = KeySwitching(self.lwe_o, self.lwe.copy(n=k * d, p=p_out, s=glwesec))
        else:
            self.lwe_o = self.lwe.copy(n=k * d, p=p_out, s=self.glwesec)

    def sample_extract(self, rotated_ct):
        N = self.d
        original_a = []
        for i in range(self.n):
            this_a = rotated_ct[i].coefficients(sparse=False)
            avec = [this_a[0]] + [-1 * this_a[N - i] for i in range(1, N)]
            original_a += avec

        b = rotated_ct[self.n].coefficients(sparse=False)[0]
        return vector(IntegerModRing(self.q), list(original_a) + [b])

    def __call__(self, c):
        """
        Switch from mod 2 to mod 3

        :param c: LWE ciphertext with plaintext modulus 2

        EXAMPLE::

            sage: msbs = ModSwitchBootstrap(LWE(n=20,q=2^10,p=2),d=64, k=3)
            sage: ct0 = msbs(msbs.lwe(0))
            sage: msbs.lwe_o.decrypt(ct0)
            0
            sage: ct1 = msbs(msbs.lwe(1))
            sage: msbs.lwe_o.decrypt(ct1)
            1

        """

        delta = self.q // self.p_out
        N = self.d
        v_modswitch = [delta] * (N // 2) + [-delta] * (N // 2)
        c = super().__call__(c, v_modswitch)
        c = self.sample_extract(c)
        if self.ksbool:
            c = self.ks(c)
            return c - vector([0] * self.lwe.n + [delta])
        else:
            return c - vector([0] * (len(c) - 1) + [delta])
