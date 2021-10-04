from sage.stats.distributions.discrete_gaussian_integer \
        import DiscreteGaussianDistributionIntegerSampler \
        as DiscreteGaussian

Zx.<x> = ZZ['x']

class RLWEDistribution:
    def __init__(self, s, N, q, sigma=3.2):
        self.n = N
        self.f = x^N + 1 # assume N = 2^k 
        self.s = s
        self.sigma = sigma
        self.D = DiscreteGaussian(sigma)
        self.Zqx = ZZ.quotient(q)['x']
        self.Rq = self.Zqx.quotient(self.f)


    def random_noise(self):
        # polinômio com grau <= n-1 e coeficientes gaussianos
        return Zx([self.D() for _ in range(self.n)])

    def random_a(self):
        return self.Rq.random_element()

    def sample(self):
        s = self.s
        a = self.random_a() # coefficientes em Zq
        e = self.random_noise()
        b = (a*s + e) # mod f e q calculado automaticamente
        return [a, b]
