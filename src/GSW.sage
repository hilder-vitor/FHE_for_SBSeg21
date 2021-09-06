load("utils.sage")
load("decompositions.sage")
load("distribution_rlwe.sage")

class GSW:
    def __init__(self, n, q, sigma, B = 2):
        self.n, self.sigma, self.B = n, sigma, B
        f = x^n + 1  # modulus defining ZZ[x] / <f>
        
        # decomposition parameters
        self.l = ceil(log(q, B))
        self.q = B^self.l
        g = Matrix(ZZ, self.l, 1, [B^i for i in range(self.l)])
        I = Matrix.identity(2)
        self.G = I.tensor_product(g) # gadget matrix

        self.Rq = ( ZZ.quotient(q) )['x'].quotient(f)
        self.sk = self.keygen()
        self.dist_rlwe = RLWEDistribution(self.sk, n, q, sigma)


    def keygen(self):
        sk = 0
        while 0 == sk:
            sk = Zx([ZZ.random_element(-1, 2) for _ in range(self.n)])
        return sk

    def encRLWE(self, m): # encrypt m into a vector ciphertext
        Rq, q, sk = self.Rq, self.q, self.sk
        a, b = self.dist_rlwe.sample()
        b += round(q / 4) * m
        return [a, b]

    def enc(self, m): # encrypt m into a matrix ciphertext
        Rq, l, sk = self.Rq, self.l, self.sk
        C = Matrix(Rq, 2*l, 2)
        for i in range(2*l):
            C[i] = self.dist_rlwe.sample()
        C += m * self.G
        return C

    def dec(self, C):
        B, l, q, sk = self.B, self.l, self.q, self.sk
        a = C[2*l - 1, 0]
        b = C[2*l - 1, 1] # == a*sk + e + (q/B) * m

        noisy_m = Zx((b - a*sk).lift())
        rounded_m = round_poly(B * noisy_m / q)
        return sym_mod_poly(rounded_m, B)

    def decRLWE(self, c):
        q, sk = self.q, self.sk
        a, b = c[0], c[1]
        noisy_m = sym_mod_poly(Zx((b - a*sk).lift()), q)
        return round_poly(4 * noisy_m / q) % 4


    ### homomorphic operations
    def add(self, C0, C1):
        C_add = C0 + C1
        return C_add

    def not_gate(self, C):
        return C + self.G

    def inv_g_row_ciphertex(self, c):
        B, l, n = self.B, self.l, self.n
        a, b = c[0], c[1]
        res = vector(self.Rq, [0] * 2 * l)
        res[0 : l] = inv_g_poly(a, B, l, n)
        res[l : 2*l] = inv_g_poly(b, B, l, n)
        return res

    def mult(self, C0, C1):
        result = Matrix(self.Rq, 2*self.l, 2)
        for i in range(2*self.l):
            decomp = self.inv_g_row_ciphertex(C0[i])
            prod_in_Rq = decomp * C1
            result[i] = prod_in_Rq
        return result

    def nand_gate(self, C0, C1):
        return self.G - self.mult(C0, C1)

    def xnor_gate(self, C0, C1):
        C00 = self.nand_gate(C0, C0)
        C11 = self.nand_gate(C1, C1)
        C0011 = self.nand_gate(C00, C11)
    
        C01 = self.nand_gate(C0, C1)
        return self.nand_gate(C0011, C01)


    def extern_prod(self, c, C):
        decomp = self.inv_g_row_ciphertex(c)
        return decomp * C
   
    ### auxiliary functions
    def get_noise_RLWE(self, c, msg):
        q, sk = self.q, self.sk
        a, b = c[0], c[1]
        b -= msg * round(q/4)
        e = b - a*sk
        e = sym_mod_poly(e, self.q)
        norm_e = infinity_norm(e)
        if norm_e == 0:
            return -1
        return log(norm_e, 2).n()

    def get_noise(self, C, msg):
        l, sk = self.l, self.sk
        C -= msg * self.G
        a = vector(C[:, 0])
        b = vector(C[:, 1]) # == a*sk + e
        e = b - a*sk
        e = sym_mod_vec(e, self.q)
        norm_e = infinity_norm_vec(e)
        if norm_e == 0:
            return -1

        return log(norm_e, 2).n()

    # defining a way to represent this class as a string
    def __str__(self):
        return "GSW: {n: %d, l: %d, B: %d, log q: %f, q: %d, sigma: %d}" % (self.n, self.l, self.B, log(self.q, 2), self.q, self.sigma)

    def size_enc_mat_MB(self):
        return (4 * self.l * self.n * self.q) / (8.0 * 10^6)


def random_poly(d, c): # return poly of degree < d and coeff. in [-c/2, c/2]
    return Zx([ZZ.random_element(-floor(c/2), ceil(c/2)) for _ in range(d)])

def test_gsw_enc_dec(gsw, num_tests=15):
    for i in range(num_tests):
        m = random_poly(gsw.n, gsw.B)
        c = gsw.enc(m)
        assert(0 == (gsw.dec(c) - m) % gsw.B)
    
    return True


def test_gsw_enc_dec_RLWE(gsw, num_tests=15):
    for i in range(num_tests):
        m = random_poly(gsw.n, 5)
        c = gsw.encRLWE(m)
        assert(0 == (gsw.decRLWE(c) - m) % 4)
    
    return True



def test_gsw_hom_add(gsw, num_tests=50):
    for i in range(num_tests):
        m1 = random_poly(gsw.n, gsw.B)
        m2 = random_poly(gsw.n, gsw.B)
        C1 = gsw.enc(m1)
        C2 = gsw.enc(m2)
        Cadd = gsw.add(C1, C2)
        assert(0 == (gsw.dec(Cadd) - (m1+m2)) % gsw.B)
    
    return True

def test_gsw_hom_mult(gsw, L=2, num_tests=10, verbose=False):
    for i in range(num_tests):
        m = random_poly(gsw.n, gsw.B)
        mand = m
        Cand = gsw.enc(m)
        for j in range(2, L+2):
            mj = random_poly(gsw.n, gsw.B)
            Cj = gsw.enc(mj)

            mand = (mand * mj) % (x^gsw.n + 1)
            Cand = gsw.mult(Cand, Cj)

            assert(0 == (gsw.dec(Cand) - mand) % gsw.B)

            if verbose:
                noisej = gsw.get_noise(Cj, mj)
                noiseand = gsw.get_noise(Cand, mand)
                tabs = "    "*(j)
                print("%snoise(c%d) = %f" % (tabs, j, noisej))
                print("%snoise(cand) = %f" % (tabs, noiseand))

        if verbose: print("")

    return True


def test_gsw_xnor(L=2, num_tests=10, verbose=False):
    sigma = 3.2
    B = 4
    l = ceil(5.5*L)
    q = B^l
    n = 4

    print("gsw = GSW(n, q, sigma, B)")
    gsw = GSW(n, q, sigma, B)
    print(gsw)


    for i in range(num_tests):
        m = ZZ.random_element(0, 2)
        mxnor = m
        Cxnor = gsw.enc(m)
        for j in range(2, L+2):
            mj = ZZ.random_element(0, 2)
            Cj = gsw.enc(mj)

            mxnor = (mxnor + mj + 1) % 2
            Cxnor = gsw.xnor_gate(Cxnor, Cj)

            assert(0 == (gsw.dec(Cxnor) - mxnor) % 2)

            if verbose:
                noisej = gsw.get_noise(Cj, mj)
                noisexnor = gsw.get_noise(Cxnor, mxnor)
                tabs = "    "*(j)
                print("%snoise(c%d) = %f" % (tabs, j, noisej))
                print("%snoise(cxnor) = %f" % (tabs, noisexnor))

        if verbose: print("")

    return True




def test_gsw_ext_prod(gsw, L=2, num_tests=10, verbose=False):
    for i in range(num_tests):
        m = random_poly(gsw.n, gsw.B)
        mand = m
        cand = gsw.encRLWE(m)
        for j in range(2, L+2):
            mj = random_poly(gsw.n, gsw.B)
            Cj = gsw.enc(mj)

            mand = (mand * mj) % (x^gsw.n + 1)
            cand = gsw.extern_prod(cand, Cj)

            assert(0 == (gsw.decRLWE(cand) - mand) % 4)

            if verbose:
                noisej = gsw.get_noise(Cj, mj)
                noiseand = gsw.get_noise_RLWE(cand, mand)
                tabs = "    "*(j)
                print("%snoise(c%d) = %f" % (tabs, j, noisej))
                print("%snoise(cand) = %f" % (tabs, noiseand))

        if verbose: print("")

    return True


def test_gsw_hom_comparison(L=3, num_tests=10, verbose=False):
    sigma = 3
    B = 2
    l = ceil(5.5*L)
    q = B^l
    n = 8

    print("gsw = GSW(n, q, sigma, B)")
    gsw = GSW(n, q, sigma, B)
    print(gsw)


    for _ in range(num_tests):
        m0 = ZZ.random_element(0, 2^L)
        m1 = ZZ.random_element(0, 2^L)
        m0 = m1 = 2^L - 1
        bits0 = m0.digits(base=2, padto=L) # lista com os n bits de m0
        bits1 = m1.digits(base=2, padto=L) # lista com os n bits de m1
        c0 = [gsw.enc(bi) for bi in bits0]
        c1 = [gsw.enc(bi) for bi in bits1]

        # compara homomorficamente
        m = 1
        c = gsw.G # enc(1) sem ruído e com a = 0
        for i in range(L):
            cmp_i = gsw.add(c0[i], c1[i]) # enc(0) <==> c0[i] == c1[i]
            cmp_i = gsw.not_gate(cmp_i) # enc(1) <==> c0[i] == c1[i]
            c = gsw.mult(c, cmp_i) # pouco ruído acumulado
            #c = gsw.mult(cmp_i, c) # muito ruído acumulado (depende de m)
            if verbose:
                m *= bits0[i] + bits1[i] + 1
                ruido_c = gsw.get_noise(c, m)
                tabs = "    "*(i)
                print("%snoise(c) = %f" % (tabs, ruido_c))
                print("%slog(m) = %f" % (tabs, log(m,2)))


        # decifra e verifica
        res = gsw.dec(c)
        assert((m0 == m1) == res)

        if verbose:
            print("")
    return True





def run_tests(lam = 50):
    n = 4
    sigma = 3
    B = 4
    l = 8
    q = B^l
    L = 3

    print("gsw = GSW(n, q, sigma, B)")
    gsw = GSW(n, q, sigma, B)

    print(gsw)

    if test_gsw_enc_dec(gsw, 25):
        print("test_gsw_enc_dec:       OK")

    if test_gsw_enc_dec_RLWE(gsw, 25):
        print("test_gsw_enc_dec_RLWE:       OK")

    if test_gsw_hom_add(gsw, 15):
        print("test_gsw_hom_add:  OK")

    if test_gsw_hom_mult(gsw, L, 5, verbose=True):
        print("test_gsw_hom_mult:  OK")

    if test_gsw_ext_prod(gsw, L, 5, verbose=True):
        print("test_gsw_ext_prod:  OK")

    L = 4
    if test_gsw_xnor(L, num_tests=4, verbose=True):
        print("test_gsw_xnor with L=%d:  OK" % L)

    L = 3
    if test_gsw_hom_comparison(L, num_tests=5, verbose=True):
        print("test_gsw_hom_comparison with %d bits:  OK" % L)


run_tests(30)
