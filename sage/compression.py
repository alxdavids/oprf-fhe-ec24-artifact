"""
FHE Ciphertext Compression.

LITERATURE:

[ACNS:CDKS21] Chen, H., Dai, W., Kim, M., & Song, Y. (2021). Efficient homomorphic
  conversion between (ring) LWE ciphertexts. In K. Sako, & N. O. Tippenhauer, ACNS 21, Part
  I (pp. 460–479). : Springer, Heidelberg.

"""
from sage.all import PolynomialRing, IntegerModRing, ZZ, vector, ceil, log, matrix
from tfhe import GLWE, KeySwitching, BlindRotation, ModSwitchBootstrap
from gadget import gadget_matrix, decompose_rlwe


class Automorphisms(GLWE):
    def __init__(self, d, q):
        """
        :param d: Degree of ring.
        :param q: Modulus.

        """
        self.d = d
        self.q = q
        self.Rq = PolynomialRing(IntegerModRing(q), "x")
        self.R = PolynomialRing(ZZ, "x")
        self.phi = self.R.gen() ** d + 1

    def __call__(self, p, t, centered=True):
        """
        Apply `X -> X^t` on polynomial ``p``.

        :param p: A polynomial or vector of polynomials
        :param t: The power to which X is raised
        :param centered: Output centered representation

        EXAMPLE::

            sage: A = Automorphisms(8,27)
            sage: p = A.Rq([1,2,3,4])
            sage: p
            4*x^3 + 3*x^2 + 2*x + 1
            sage: A(p,3)
            3*x^6 + 2*x^3 - 4*x + 1

            sage: q = A.Rq([0]*4+[1,2,3,4])
            sage: A([p,q],3)
            (3*x^6 + 2*x^3 - 4*x + 1, -2*x^7 + 4*x^5 - x^4 + 3*x^2)
        """
        try:
            _ = len(p)
            return vector([self.__call__(pi, t, centered) for pi in p])
        except TypeError:
            r = self.Rq((p(self.R.gen() ** t) % self.phi))
            if centered:
                return self.lift_centered(r)
            else:
                return r


class AutoKeySwitching(KeySwitching):
    def __init__(self, rlwe, B=2):
        self.B = B
        self.ell = ceil(log(rlwe.q, B))
        self.rlwe = rlwe
        self.autos = Automorphisms(rlwe.d, rlwe.q)
        self.auto_s = []
        self.auto_ksk = self.key_gen()

    def key_gen(self):
        ak = []
        rlwe = self.rlwe

        self.lwe_o = rlwe
        powers = [rlwe.d / 2**t + 1 for t in range(0, log(rlwe.d, 2))]
        for power in powers:
            self.lwe_i = GLWE(
                d=rlwe.d, q=rlwe.q, k=rlwe.n, p=rlwe.p, s=[self.autos(rlwe._s[0], power)]
            )
            akt = super().key_gen()[0]
            self.auto_s.append(self.autos(rlwe._s[0], power))
            ak.append(akt)
        return ak

    def eval_auto(self, c, t):
        c_auto = self.autos(c, t, False)
        a, b = self.rlwe.Rq(c_auto[0]), self.rlwe.Rq(c_auto[1])

        ctxt_out = vector([0] * self.rlwe.n + [b])
        adec = decompose_rlwe([a], self.B, self.ell, self.rlwe.d)
        t_ind = ZZ(log(self.rlwe.d / (t - 1), 2))
        for j in range(self.ell):
            ctxt_out += (adec[j] * self.auto_ksk[t_ind][j]) % self.rlwe.phi

        return ctxt_out

    def trace_eval(self, c, n=1):
        """
        Evaluate homomorphic trace i.e. zero all non-constant coeffs (adding a factor d)

        :param c: A rlwe ciphertext
        :param n: If n>1 only clear part of non-constant coefficients with smaller factor

        EXAMPLE::

            sage: glwe = GLWE(4, 2**32-1, 1, p=16, s=[[1,1]])
            sage: c0 = glwe([7, 1,1,1])
            sage: AK = AutoKeySwitching(glwe)
            sage: glwe.decrypt(AK.trace_eval(c0))
            12
            sage: (7 * 4) % 16
            12

            sage: glwe.decrypt(AK.trace_eval(c0, 2))
            2*x^2 + 14

        """
        ctxt = c
        rlwe = self.rlwe
        powers = [rlwe.d / 2**t + 1 for t in range(0, log(rlwe.d / n, 2))]
        for power in powers:
            c_auto = self.eval_auto(ctxt, power)
            ctxt += c_auto

        return vector(self.rlwe.Rq, ctxt)

    def coeff_extract(self, c):
        """
        Extract d rlwe ciphertexts from a single rlwe ciphertext. The d ciphertexts
        encrypt the coefficients of the input rlwe ciphertext with a factor of d.

        :param c: A rlwe ciphertext to extract

        EXAMPLE::

            sage: glwe = GLWE(4, 2**32-1, 1, p=16)
            sage: c0 = glwe([2,1,3,4])
            sage: AK = AutoKeySwitching(glwe)
            sage: ct_list = AK.coeff_extract(c0)
            sage: [glwe.decrypt(ct_list[i]) for i in range(4)]
            [8, 4, 12, 0]


        """
        rlwe = self.rlwe
        res = [vector([]) for _ in range(rlwe.d)]
        res[0] = c
        powers = [rlwe.d / 2**t + 1 for t in range(log(rlwe.d, 2))]
        for power in powers:
            for j in range(rlwe.d / (power - 1)):
                old = res[j]
                tmp = self.eval_auto(old, power)
                res[j] = old + tmp
                poX = -1 * self.rlwe.Rq.gen() ** (rlwe.d - rlwe.d / (power - 1)) % self.rlwe.phi
                res[j + rlwe.d / (power - 1)] = poX * (old - tmp) % self.rlwe.phi

        return res


class PackedBlindRotation(BlindRotation):
    """
    Blind rotation via a compressed blind rotation key. The constructor generates
    a compressed key and then decompresses it.

    EXAMPLE::

        sage: from tfhe import LWE
        sage: pbr = PackedBlindRotation(LWE(n=2, q=2^15-1, p=4), d=32)
        sage: c1 = pbr(pbr.lwe(1))
        sage: GLWE.decrypt(pbr, c1)[0]
        1
        sage: c0 = pbr(pbr.lwe(0))
        sage: GLWE.decrypt(pbr, c0)[0]
        0
        sage: all([GLWE.decrypt(pbr, pbr(pbr.lwe(1)))[0] == 1 for _ in range(32)])
        True
        sage: all([GLWE.decrypt(pbr, pbr(pbr.lwe(0)))[0] == 0 for _ in range(32)])
        True

    """

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
        :param s: Explicitly set a" secret.

        """
        super(BlindRotation, self).__init__(d, lwe.q, k, chi_s, chi_e, B, p, s)
        self.lwe = lwe
        self.rlwe = self.base_scheme()
        self.auto_ks = AutoKeySwitching(self.rlwe, B)
        self.sqk = self.create_square_key()
        self.compressed_brk = self.key_gen()
        self.brk = self.decompress()

    def key_gen(self):
        full_list = []
        lwe_n = self.lwe.n
        invN = self.Rq(1 / (self.d))
        gadget_vec = gadget_matrix(1, self.B, self.ell)
        for i in range(lwe_n):
            full_list += list(invN * self.lwe._s[i] * gadget_vec)[0]

        brk_list = []
        for i in range((self.ell * lwe_n) // self.d + 1):
            this_poly = self.rlwe.Rq(full_list[i * self.d : (i + 1) * self.d])
            brk_list.append(self.rlwe(this_poly, raw=True))

        return brk_list

    def create_square_key(self):
        rlwe = self.rlwe
        square = (rlwe._s[0]) ** 2 % self.phi
        square_key = []
        for i in range(self.ell):
            square_key.append(rlwe(square * self.B**i, raw=True))

        return matrix(square_key)

    def eval_square_mult(self, ct):
        """
        Send a ciphertext encrypting `m` to a ciphertext encrypting `s⋅m`

        :param ct: RLWE ciphertext

        EXAMPLE::

            sage: from tfhe import LWE
            sage: lwe = LWE(4, 2**15-3, p=16)
            sage: pbr = PackedBlindRotation(lwe, 4, s=[(1,0,1,1)], p=16)
            sage: pbr.rlwe._s[0]
            x^3 + x^2 + 1
            sage: c0 = pbr.rlwe([0,2,0,3])
            sage: pbr.rlwe.decrypt(pbr.eval_square_mult(c0))
            5*x^3 + 13*x^2 + 15*x + 14
            sage: (pbr.rlwe.decrypt(c0) * pbr.rlwe._s[0]) % pbr.rlwe.phi
            5*x^3 - 3*x^2 - x - 2

        """
        a, b = ct[0], ct[1]
        adec = decompose_rlwe([a], self.B, self.ell, self.d)
        return vector((vector([b, 0]) + adec * self.sqk) % self.rlwe.phi)

    def decompress(self):
        cbrk_list = self.compressed_brk
        ct_list = []
        for ct in cbrk_list:
            ct_list += self.auto_ks.coeff_extract(ct)
        # encryption of si*B**j is i*ell+jth element
        sqct_list = []
        for ct in ct_list:
            sqct_list.append(self.eval_square_mult(ct))

        # now arrange as RGSW ctxts of si
        brk_list = []
        ell = self.ell
        for i in range(self.lwe.n):
            rgswi = sqct_list[i * ell : (i + 1) * ell] + ct_list[i * ell : (i + 1) * ell]
            brk_list.append(matrix(rgswi))

        return tuple(brk_list)


class PackedModSwitchBootstrap(ModSwitchBootstrap, PackedBlindRotation):
    pass


class CiphertextCompression(AutoKeySwitching):
    def __init__(self, lwe, d, B=2):
        self.lwe = lwe
        self.autos = Automorphisms(d, lwe.q)
        self.rlwe = GLWE(d, lwe.q, p=lwe.p, s=[self.convert_secret(d, lwe._s[: lwe.n])])
        super().__init__(self.rlwe, B)

    def convert_secret(self, d, lwe_s):
        """
        Convert an LWE secret to a RLWE secret for the LWE to RLWE conversion.
        This conversion means the same automorphism key is used for key and
        ciphertext compression.

        :param d: The ring dimension
        :param lwe_s: An lwe secret

        EXAMPLE::

            sage: from tfhe import LWE
            sage: lwe = LWE(4, 2**15-3, p=16, s=(1,0,1,0))
            sage: cc = CiphertextCompression(lwe, 4)
            sage: cc.convert_secret(4,lwe._s[:lwe.n])
            x^2 + 1
            sage: lwe._s[:lwe.n]
            (1, 0, 1, 0)

        """
        R = self.autos.R
        ring_s = R(list(lwe_s))
        return ring_s

    def pack_lwes(self, ct_list):
        """
        Pack a list of rlwe ciphertexts (output from to_rlwe function) into a single
        rlwe ciphertext.

        EXAMPLE::

            sage: from tfhe import LWE
            sage: lwe = LWE(4, 2**15-3, p=16)
            sage: cc = CiphertextCompression(lwe, 4)
            sage: ct_list = [cc.to_rlwe(lwe(i)) for i in range(4)]
            sage: cc.rlwe.decrypt(cc.pack_lwes(ct_list))
            12*x^3 + 8*x^2 + 4*x

        """
        ell = log(len(ct_list), 2)
        N = self.rlwe.d
        if ell == 0:
            return ct_list[0]
        else:
            ct_even = self.pack_lwes([ct_list[2 * j] for j in range(2 ** (ell - 1))])
            ct_odd = self.pack_lwes([ct_list[2 * j + 1] for j in range(2 ** (ell - 1))])
            poX = self.rlwe.Rq.gen() ** (N / (2**ell)) % self.rlwe.phi
            term1 = (ct_even + poX * ct_odd) % self.rlwe.phi
            term2 = self.eval_auto((ct_even - poX * ct_odd) % self.rlwe.phi, 2**ell + 1)
            ct = (term1 + term2) % self.rlwe.phi
            return ct

    def to_rlwe(self, lwe_ct):
        """
        Embed a lwe ciphertext into an rlwe one. The output ciphertext encrypts the
        lwe plaintext in its constant. This uses a LWE to RLWE TFHE transform
        so that the same automorphism key can be used for key compression and
        ciphertext compression.

        EXAMPLE::

            sage: from tfhe import LWE
            sage: lwe = LWE(4, 2**15-3, p=16)
            sage: cc = CiphertextCompression(lwe, 4)
            sage: cc.rlwe.decrypt(cc.to_rlwe(lwe(3)))[0]
            3

        """
        n = self.lwe.n
        d = self.rlwe.d
        a, b = list(lwe_ct[:n]) + [0] * (d - n), lwe_ct[n]
        R = self.autos.R
        phi = R.gen() ** d + 1
        ring_a = R(list(a))
        a = R((ring_a(-R.gen() ** (d - 1)) % phi))
        return vector([self.rlwe.Rq(a), self.rlwe.Rq(b)])

    def lwes_to_rlwe(self, ct_list):
        """
        Use ``to_rlwe`` and ``pack_lwes`` to pack many LWE ciphertexts into a single RLWE
        ciphertext. The output plaintext contains the lwe plaintexts in the non-zero
        coefficients (multiplied by N mod p).

        :param ct_list:

        EXAMPLE::

            sage: from tfhe import LWE
            sage: lwe = LWE(8, 2**15, p=3)
            sage: cc = CiphertextCompression(lwe, 8)
            sage: cc.rlwe.decrypt(cc.lwes_to_rlwe(list(map(lwe,[1,1,1,1, 1,1,1,1]))))
            2*x^7 + 2*x^6 + 2*x^5 + 2*x^4 + 2*x^3 + 2*x^2 + 2*x + 2
            sage: cc.rlwe.decrypt(cc.lwes_to_rlwe(list(map(lwe,[0,1,1,1, 1,1,1,1]))))
            2*x^7 + 2*x^6 + 2*x^5 + 2*x^4 + 2*x^3 + 2*x^2 + 2*x
            sage: cc.rlwe.decrypt(cc.lwes_to_rlwe(list(map(lwe,[0,1,0,1, 0,1,0,1]))))
            2*x^7 + 2*x^5 + 2*x^3 + 2*x
            sage: cc.rlwe.decrypt(cc.lwes_to_rlwe(list(map(lwe,[0,2,0,2, 0,2,0,2]))))
            x^7 + x^5 + x^3 + x

        """
        n = len(ct_list)
        rlwe_ct_list = []
        for ct in ct_list:
            rlwe_ct_list.append(self.to_rlwe(ct))

        ct = self.pack_lwes(rlwe_ct_list)

        return self.trace_eval(ct, n)
