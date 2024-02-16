"""
OPRF Candidate from TFHE and Crypto Dark Matter PRF Candidate.

Implements:

- OPRF
- circuit privacy
- ciphertext compression
- bootrstrapping key compression

Does not implement:

- VOPRF/check points
- NIZK proofs

"""

from sage.all import (
    GF,
    codes,
    gcd,
    identity_matrix,
    random_matrix,
    set_random_seed,
    matrix,
    vector,
    ZZ,
)
from tfhe import ModSwitchBootstrap, apply_plaintext_matrix
from cpbs import CPModSwitchBootstrap, CPPackedModSwitchBootstrap
from compression import PackedModSwitchBootstrap, CiphertextCompression


class WeakPRF:
    """
    Based on Construction 3.1 of:

    - Boneh, D., Ishai, Y., Passelègue, A., Sahai, A., & Wu, D. J. (2018).
      Exploring crypto dark matter: new simple PRF candidates and their
      applications. Cryptology ePrint Archive, Report 2018/1218.
      https://eprint.iacr.org/2018/1218
    """

    def __init__(self, m_p=256, n_p=256, m_bound=128, p=2, q=3, t=None, seed=None):
        # p.40, optimistic, λ=128
        self.m_p, self.n_p = m_p, n_p
        if seed is not None or t is not None:
            set_random_seed(hash((t, seed)))
        self.A = random_matrix(GF(p), m_p, n_p + 1)
        self.q = q

        if gcd(p, q) != 1:
            raise ValueError(f"p={p} and q={q} must be coprime.")

        if m_bound == 1:
            self.G_out = matrix(GF(q), 1, m_p, [1] * m_p)
        else:
            i = 2
            for i in range(2, m_p // 2):
                C = codes.BCHCode(GF(q), m_p, i)
                if C.dimension() <= m_bound:
                    self.G_out = C.generator_matrix()
                    break
            else:
                raise RuntimeError

    def __call__(self, x):
        y = self.A * vector(list(x) + [1])
        y = y.lift_centered()
        z = self.G_out * y
        return z


class PRF:
    def __init__(self, m_p=256, n_p=256, n=128, m_bound=128, p=2, q=3, t=None, seed=None):
        if p != 2 or q != 3:
            raise NotImplementedError
        self.F = WeakPRF(m_p=m_p, n_p=n_p, m_bound=m_bound, p=p, q=q, t=t, seed=seed)
        self.n = n
        self.G_inp = identity_matrix(GF(q), n).stack(
            random_matrix(GF(q), (n_p - n) // 2, n)
        )

    def __call__(self, x):
        x = vector(x).lift().change_ring(GF(self.F.q))
        y = self.G_inp * x
        z = []
        for y_ in y[: self.n].lift():
            z.append(y_)
        for y_ in y[self.n :].lift():
            z.append(y_ % 2)
            z.append(y_ // 2)
        z = vector(GF(2), self.F.n_p, z)
        return self.F(z)


class OPRF(PRF):
    """
    EXAMPLE::

        sage: set_random_seed(1337)
        sage: from tfhe import LWE
        sage: from oprf import OPRF
        sage: oprf = OPRF(LWE(4, 3*7681, "binary", 3.0, p=2))
        sage: oprf([0]*8)
        (1, 0, 2, 2, 1, 1, 0, 2, 2, 1, 2, 2)
        sage: c = oprf.blind_eval([0]*8)
        sage: vector([oprf.msbs.lwe_o.decrypt(c_) for c_ in c])
        (1, 0, 2, 2, 1, 1, 0, 2, 2, 1, 2, 2)
    """

    def __init__(
        self,
        lwe,
        d=256,
        k=1,
        m_p=16,
        n_p=16,
        n=8,
        m_bound=16,
        p=2,
        q=3,
        t=None,
        seed=None,
        cp=False,
        key_pack=False,
        ct_pack=False,
        B=128,
    ):
        """

        :param lwe: LWE instance for blind evaluations
        :param d: MLWE ring dimension.
        :param k: MLWE module rank.
        :param m_p: Number of rows of `A`.
        :param n_p: Number of columns of `A`.
        :param n: Input dimension of vector mod `p`
        :param m_bound: Bound on the output dimension of vector mod `q`
        :param p: Input modulus
        :param q: Output modulus ≠ p
        :param t: Plain input
        :param seed: randomness seed
        :param cp: Enable circuit privacy (extremely slow)
        :param key_pack: Public-key compression
        :param ct_pack: RLWE packing
        :param B: Gadget decomposition base.

        """
        super().__init__(m_p=m_p, n_p=n_p, n=n, m_bound=m_bound, p=p, q=q, t=t, seed=seed)
        self.lwe = lwe

        if not ZZ(3).divides(lwe.q):
            raise ValueError(f"3∤{lwe.q}.")

        self.d = d
        if cp and not key_pack:
            self.msbs = CPModSwitchBootstrap(self.lwe, d=d, k=k, B=B)
        elif cp and key_pack:
            self.msbs = CPPackedModSwitchBootstrap(self.lwe, d=d, k=k, B=B)
        elif (not cp) and key_pack:
            self.msbs = PackedModSwitchBootstrap(self.lwe, d=d, k=k, B=B)
        else:
            self.msbs = ModSwitchBootstrap(self.lwe, d=d, k=k, B=B)

        self.ct_pack = ct_pack
        if ct_pack:
            self.packing = CiphertextCompression(self.msbs.lwe_o.copy(p=q), self.msbs.d)

    def blind_eval(self, x):
        n = self.n
        x = vector(x).lift().change_ring(GF(self.F.q))
        y = self.G_inp * x
        z = []
        for y_ in y[:n].lift():
            z.append(y_)
        for y_ in y[n:].lift():
            z.append(y_ % 2)
            z.append(y_ // 2)
        z = vector(GF(2), self.F.n_p, z)

        c = [self.lwe(z_) for z_ in z] + [self.lwe(1)]

        # we might as well randomize
        c = apply_plaintext_matrix(self.F.A, c, randomize=True)
        self.prebs = c

        c = [self.msbs(c_) for c_ in c]
        c = apply_plaintext_matrix(self.F.G_out, c)

        if self.ct_pack:
            packed = self.packing.lwes_to_rlwe(c)
            return packed
        else:
            return c
