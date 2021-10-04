load("utils.sage")
load("decompositions.sage")
load("distribution_rlwe.sage")

class GSW:
    def __init__(self, n, q, sigma, B = 2):
        self.n, self.sigma, self.B = n, sigma, B
        f = x^n + 1  # modulus que define ZZ[x] / <f>
        
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

    def encRLWE(self, m): # encrypt m into a matrix ciphertext
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
