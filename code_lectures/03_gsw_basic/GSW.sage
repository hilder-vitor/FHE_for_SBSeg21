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

    ### homomorphic operations
#    def add(self, C0, C1):
        # implementar durante o curso

    # defining a way to represent this class as a string
    def __str__(self):
        return "GSW: {n: %d, l: %d, B: %d, log q: %f, q: %d, sigma: %d}" % (self.n, self.l, self.B, log(self.q, 2), self.q, self.sigma)
