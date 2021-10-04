from sage.stats.distributions.discrete_gaussian_integer \
        import DiscreteGaussianDistributionIntegerSampler \
        as DiscreteGaussian

class LWEDistribution:
    def __init__(self, s, q, sigma=3.2):
        self.n = len(s)
        self.s = s
        self.sigma = sigma
        self.D = DiscreteGaussian(sigma)
        self.Zq = ZZ.quotient(q)

    def random_noise(self):
        return self.D()

    def random_a(self):
        a = [self.Zq.random_element() for _ in range(self.n)]
        a = vector(self.Zq, a) # converte lista em vetor
        return a

    def sample(self):
        n, s = self.n, self.s
        a = self.random_a() # vetor em Zq^n
        e = self.random_noise()
        b = a*s + e
        return a, b
