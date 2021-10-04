load("lwe_base_scheme.sage")
load("GSW.sage")

def is_binary(v):
    for vi in v:
        if vi != 0 and vi != 1:
            return False
    return True


def round_vec(v):
    return vector([round(vi) for vi in v])

#   Given a polynomial g, this function returns a vector representing g % f.
#   The dimension of the returned polynomial is deg(f).
#   Notice that this map is an isomorphism from ZZ[x]/<x^N+1> to ZZ^N.
def poly_to_vec(g, N):
    f = x^N+1
    v = (g % f).coefficients(sparse=False)
    return vector(v + [0]*(N - len(v)))

#   Given a polynomial g, this function returns a matrix that represents
# the product by g in the ring ZZ[x]/<x^N+1>, that is,
# a matrix A that satisfies A * poly_to_vec(h) = poly_to_vec(h*g % f) for all h.
#   The returned matrix is NxN.
def poly_to_mat(g, N):
    A = Matrix(ZZ, [poly_to_vec(g * x^i, N) for i in range(N)])
    #return A.transpose()
    return A


class Bootstrapper:
    # _gsw scheme used as the accumulator
    # lwe_base_scheme is the scheme used to generate the ciphertexts we want to
    #    refresh
    def __init__(self, _gsw, lwe_scheme):
        assert(_gsw.q == lwe_scheme.q)
        self.gsw = _gsw
        self.base_scheme = lwe_scheme
        self.one = self.gsw.G
        self.boot_key_gen()
        self.key_swt_gen()

    # sk: secret key of the base scheme
    def boot_key_gen(self):
        s = self.base_scheme.sk
        assert(is_binary(s))
        n = len(s)
        self.bk = [self.gsw.enc(s[i]) for i in range(n)]

    def key_swt_gen(self):
        N, Q = self.gsw.n, self.gsw.q
        q = self.base_scheme.q
        l = ceil(log(Q, 2))
        z = poly_to_vec(self.gsw.sk, N)

        # K_{i, j} = enc(s_i * 2^j)
        K = [0] * N
        for i in range(N):
            zi = z[i]
            enc_zi = [0] * l
            for j in range(l):
                c = self.base_scheme.enc(0)
                c[1] = (c[1] + 2^j * zi) % q
                enc_zi[j] = c
            K[i] = enc_zi
        self.ksk = K

    # Receives an LWE ciphertext (a, b) encrypted under the key
    # poly_to_vec(gsw.sk) and outputs one encrypted under the
    # key of the base scheme
    def switch_key(self, c):
        a, b = c
        q = self.base_scheme.q
        l = ceil(log(q, 2))
        n = self.base_scheme.n
        K = self.ksk
        out_a = vector([0]*n)
        out_b = 0
        for i in range(len(a)):
            v = vector(ZZ(a[i]).digits(base=2, padto=l))
            for j in range(l):
                if 1 == v[j]:
                    out_a += K[i][j][0] # adding "a" term of key
                    out_b += K[i][j][1] # adding "b" term of key

        # Now out_b = out_a * s + e + a * z
        # so, b - out_b = -out_a*s - e + e_b + (q / 4) * m
        out_b = (b - out_b) % q
        out_a = (-out_a ) % q
        return [out_a, out_b]


    #   Receives an RLWE ciphertext acc encrypting some message m, 
    # an integer a_i with absolute value smaller than or equal to 2*N, 
    # and a GSW ciphertext bk_i encrypting a bit s_i.
    #   Returns an RLWE ciphertext encrypting m * x^(-a_i * s_i)
    def cmux_gate(self, acc, a_i, bk_i):
        N = self.gsw.n
        C = self.one + (x^(2*N-a_i) - 1) * bk_i
        acc = self.gsw.extern_prod(acc, C)
        return acc

    # Receives an LWE ciphertext (a, b) in Z_q^(n+1) and returns an LWE 
    # ciphertext defined modulo 2*N encrypting the same message
    def mod_switch_to_2N(self, a, b):
        q, N = self.base_scheme.q, self.gsw.n
        _a = vector([round(ZZ(ai) * 2*N/q) for ai in a])
        _b = round(ZZ(b) * 2*N / q)
        return _a, _b

    # receives 
    #       a base-scheme ciphertext c = (a, b) in Z_q^(n+1)
    #       the bootstrapping key bk
    # and homomorphically computes the linear part of the decryption,
    # i.e., b - a*sk, generating an RLWE encryption of X^(e + Delta*m)
    def linear_part_dec(self, c):
        Q, N, RQ = self.gsw.q, self.gsw.n, self.gsw.Rq
        a, b = c
        acc = [RQ(0), round(Q/8) * RQ(x)^(N // 2 + b)] # trivial, noiseless, RLWE enc of x^b
        for i in range(len(a)):
            acc = self.cmux_gate(acc, a[i], self.bk[i])
        return acc

    def convert_rlwe_lwe(self, acc):
        N, Q = self.gsw.n, self.gsw.q
        a, b = Zx(acc[0].lift()), Zx(acc[1].lift())

        u = vector([-1] * N)

        lwe_a = (poly_to_mat(a, N) * u) % Q
        lwe_b = poly_to_vec(b, N) * u # == lwe_a * coefs(s) + e + (Q/8)*(2*m - 1)
        lwe_b = (lwe_b + round(Q/8)) % Q # == lwe_a * coefs(s) + e + (Q/4) * m
        
        return lwe_a, lwe_b # in Z_Q^(N+1)

    def refresh(self, c):
        a, b = self.mod_switch_to_2N(c[0], c[1])
        acc = self.linear_part_dec([a, b])

        c = self.convert_rlwe_lwe(acc)
    
        c = self.switch_key(c)

        return c


def test_bootstrapp(level=4):
    n, N, sigma, = 2^3, 2^5, 3.2
    B, l = 2, 17
    q = Q = B^l # assume q = Q

    lwe_scheme = LWEScheme(n, q, sigma)
    print(lwe_scheme)
    gsw = GSW(N, Q, sigma, B)
    print(gsw)

    boot = Bootstrapper(gsw, lwe_scheme)

    m = random_bit()
    c = lwe_scheme.enc(m)

    for i in range(1, level+1):
        mi = random_bit()
        m = (1 - m * mi)
        ci = lwe_scheme.enc(mi)
        cnand = lwe_scheme.nand(c, ci)

        print("nand(m, m%d) = %s" % (i, lwe_scheme.dec(cnand)))

        c = boot.refresh(cnand)

        dec_m = lwe_scheme.dec(c, level=1)
    
        assert(m == dec_m)

        noise = lwe_scheme.get_noise(c, m, level=1)
        print("Refreshed noise: %f" % noise)


test_bootstrapp()
