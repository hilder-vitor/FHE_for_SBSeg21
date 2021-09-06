load("utils.sage")
load("distribution_acd.sage")

#    A simple implementation of the ACD-based homomorphic encryption scheme 
# presented in the paper ''Fully Homomoprhic Encryption over the Integers'', 
# by Dijk,  Gentry,  Halevi,  and  Vaikuntanathan, in Eurocrypt 2010.
#    This implementation does not include public keys nor bootstrapping.
class DGHV:

    def __init__(self, gamma, eta, rho, t = 2, p = 1):
        assert(gamma > eta)
        assert(eta > rho)
        if 1 == p:
            p = random_prime(2^eta, lbound=2^(eta - 1))
        else:
            assert(eta-1 <= p.nbits() <= eta) # if p is given, it must have eta bits

        self.gamma = gamma
        self.eta = eta
        self.rho = rho
        self.t = t
        self.p = p
        self.x0 = p * sample_q(gamma, eta) #+ self.sample_r()
        self.Zp = ZZ.quotient(p)
        self.Zx0 = ZZ.quotient(self.x0)
        
    def enc(self, m):
        q = sample_q(self.gamma, self.eta)
        r = sample_r(self.rho)
        c = self.p*q + self.t * r + m
        c %= self.x0
        return c

    # c is a polynomial in R = Z[x] / <x^N - 1>
    def dec(self, c):
        noisy_msg = sym_mod(c, self.p) # == t * r + msg
        return noisy_msg % self.t

    ### homomorphic operations
    def not_gate(self, c):
        return (1 - c) #% self.x0

    def add(self, c1, c2):
        return (c1 + c2) #% self.x0

    def mult(self, c1, c2):
        return c1 * c2 #%self.x0


    ### auxiliary functions
    def get_noise(self, c, msg):
        noisy_msg = sym_mod(c, self.p) # == t * r + msg
        noise = abs((noisy_msg - msg) / self.t)
        if (noise > 0):
            return log(noise, 2).n()
        else:
            return 0
    
    # defining a way to represent this class as a string
    def __str__(self):
        return "DGHV: {gamma: %d, eta: %d, rho: %d, t: %d,  p: %d}" \
                 % (self.gamma, self.eta, self.rho, self.t, self.p)


def random_bit():
    if random() > 0.5:
        return 1
    return 0

def test_dghv_enc_dec(dghv, num_tests=15):
    for i in range(num_tests):
        m = ZZ.random_element(0, dghv.t)
        c = dghv.enc(m)
        assert(dghv.dec(c) == m)
    
    return True

def test_dghv_hom_not_gate(dghv, num_tests=50):
    for i in range(num_tests):
        m = random_bit()
        c = dghv.enc(m)
        notc = dghv.not_gate(c)
        assert(dghv.dec(notc) == (not m))
    return True


def test_dghv_hom_add(dghv, num_tests=50):
    t = dghv.t
    for i in range(num_tests):
        m1 = ZZ.random_element(0, t)
        m2 = ZZ.random_element(0, t)
        c1 = dghv.enc(m1)
        c2 = dghv.enc(m2)
        cadd = dghv.add(c1, c2)
        assert(0 == (dghv.dec(cadd) - (m1+m2)) % dghv.t)
    
    return True

def test_dghv_hom_mult(dghv, L=2, num_tests=10, verbose=False):
    rho, eta, gamma, t = dghv.rho, dghv.eta, dghv.gamma, dghv.t
    for i in range(num_tests):
        m = ZZ.random_element(0, t)
        mand = m
        cand = dghv.enc(m)
        for j in range(2, L+2):
            mj = ZZ.random_element(0, t)
            cj = dghv.enc(mj)

            mand = (mand * mj) % t
            cand = dghv.mult(cand, cj)

            assert(dghv.dec(cand) == (mand))

            if verbose:
                noisej = dghv.get_noise(cj, mj)
                noiseand = dghv.get_noise(cand, mand)
                tabs = "    "*(j)
                print("%snoise(c%d) = %f" % (tabs, j, noisej))
                print("%snoise(cand) = %f" % (tabs, noiseand))
                print("%sbit length cand = %d" % (tabs, log(abs(cand),2)))
                print("%sgamma = %d" % (tabs, dghv.gamma))

        if verbose: print("")

    return True


def test_dghv_hom_comparison(dghv, n=3, num_tests=10, verbose=False):
    for _ in range(num_tests):
        m0 = ZZ.random_element(0, 2^n)
        m1 = ZZ.random_element(0, 2^n)
        bits0 = m0.digits(base=2, padto=n) # lista com os n bits de m0
        bits1 = m1.digits(base=2, padto=n) # lista com os n bits de m1
        c0 = [dghv.enc(bi) for bi in bits0]
        c1 = [dghv.enc(bi) for bi in bits1]

        # compara homomorficamente
        c, m = 1, 1
        for i in range(n):
            cmp_i = dghv.add(c0[i], c1[i]) # enc(0) <==> c0[i] == c1[i]
            cmp_i = dghv.not_gate(cmp_i) # enc(1) <==> c0[i] == c1[i]
            c = dghv.mult(c, cmp_i) # c *= cmp_i
            if verbose:
                m *= ((bits0[i] + bits1[i] + 1) % 2)
                ruido_c = dghv.get_noise(c, m)
                tabs = "    "*(i)
                print("%snoise(c) = %f" % (tabs, ruido_c))
                print("%sbit length c = %d" % (tabs, log(abs(c),2)))
                print("%sgamma = %d" % (tabs, dghv.gamma))


        # decifra e verifica
        res = dghv.dec(c)
        assert((m0 == m1) == res)

        if verbose:
            print("")
    return True




def run_tests(lam = 50):
    rho = lam
    print("rho = %d" % rho)
    L = 4
    t = 2
    eta = (L+3)*lam
    print("eta = %d" % eta)
    gamma = max(ceil(lam * (eta - rho)^2 / (log(lam, 2))), 2*eta)
    print("gamma = %d" % gamma)

    dghv = DGHV(gamma, eta, rho, t)
    print(dghv)

    if test_dghv_enc_dec(dghv, 25):
        print("test_dghv_enc_dec:       OK")

    if test_dghv_hom_not_gate(dghv, 15):
        print("test_dghv_hom_not_gate:  OK")

    if test_dghv_hom_add(dghv, 15):
        print("test_dghv_hom_add:  OK")

    if test_dghv_hom_mult(dghv, L, 5, verbose=True):
        print("test_dghv_hom_and_gate:  OK")

    if test_dghv_hom_comparison(dghv, n=L, num_tests=20, verbose=True):
        print("test_dghv_hom_comparison:  OK")

run_tests(30)

