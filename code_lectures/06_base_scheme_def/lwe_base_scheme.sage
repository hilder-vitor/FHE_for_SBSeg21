load("utils.sage")
load("distribution_lwe.sage")

class LWEScheme:
    def __init__(self, n, q, sigma=3.2):
        self.n, self.q = n, q
        self.sigma = sigma
        self.D = DiscreteGaussian(sigma)
        self.Zq = ZZ.quotient(q)
        self.keygen()
        self.dist_lwe = LWEDistribution(self.sk, q, sigma)


    def keygen(self):
        n = self.n
        bin_list = [ZZ.random_element(0, 2) for _ in range(n)]
        self.s = vector(ZZ, bin_list)
        self.sk = self.s

    def enc(self, m):
        a, _b = self.dist_lwe.sample()
        b = _b + round(self.q / 4) * m
        return [a, b]

    # level 1: q/4, erro(c) < q/16
    def dec(self, c):
        q, s = self.q, self.s
        a, b = c
        # then b = a*s + e + (q/4)*m and abs(e) < q/16
        noise_bound = round(q/16)
        Delta = round(q/4)
        noisy_m = (b - a*s + noise_bound) # == e' + Delta*m with |e'| < Delta
        if noisy_m < Delta:
            return 0
        else:
            return 1

    def get_noise(self, c, m, level):
        q = self.q
        Delta = q // 4
        noise = sym_mod(c[1] - c[0]*self.sk - Delta*m, q)
        if 0 == abs(noise):
            return 0
        return log(abs(noise), 2)

    def __str__(self):
        return "LWE scheme: {n: %d, log q: %f," \
            "sigma: %f, q: %d}" \
         % (self.n, log(self.q,2), \
         self.dist_lwe.sigma, self.q)
