#################################################################################
#       Orthogonal lattice attack on the original AGCD
#   as described in the paper Algorithms for the approximate common
#   divisor problem, by Galbraith, Gebregiyorgis and Murphy.
#
#   No noiseless x0
#   Single prime
#
#       Notice that you must by yourself the block size of BKZ that achieves the 
#   chosen root-Hermite factor.
#################################################################################

from sage.modules.free_module_integer import IntegerLattice as Lattice


def sample_r(rho):
	return ZZ.random_element(-2^rho+1,2^rho)

def sample_vec_r(rho, m):
	return vector(ZZ, [sample_r(rho) for _ in range(m)])

def sample_q(gamma, eta):
	return ZZ.random_element(0, 2^(gamma - eta))

def sample_vec_q(gamma, eta, m):
	return vector(ZZ, [sample_q(gamma, eta) for _ in range(m)])

def sample_p(eta):
    return random_prime(2^eta, False, 2^(eta-1))

def sample_agcd(gamma, p, rho):
    eta = ceil(log(p, 2))
    r = sample_r(rho)
    q = sample_q(gamma, eta)
    c = p*q + r
    return c



def sym_mod(a, n):
    a = a % n
    if 2*a > n:
        return a - n
    return a



def orthogonal_lattice_atk():

    # Function that uses t AGCD samples to generate
    # a matrix used as the lattice basis
    def sample_basis(gamma, p, rho, t):
        R = 2^rho
        B = Matrix(ZZ, t, t+1)
        for i in range(t):
            B[i, 0] = sample_agcd(gamma, p, rho)
            B[i, i+1] = R
        return B 

    rho = 40
    eta = 50
    rhfLLL = 1.014
    gamma_discriminant_neg = ZZ(ceil((eta - rho)^2 / (4 * log(rhfLLL,2))));
    gamma = ceil(gamma_discriminant_neg / 2)

    print("gamma to rule out LLL: %d" % gamma_discriminant_neg)
    print("chosen gamma: %d" % gamma)

    t = ceil(sqrt( gamma / (4*log(rhfLLL,2)) ))

    p = sample_p(eta)

    B = sample_basis(gamma, p, rho, t)


    L = Lattice(B)
    L.LLL()
#    L.BKZ(block_size=30, prune=15) # XXX: you must choose the correct block size for the given root-Hermite factor


    total = 0
    vecs_ort_q = []
    avg_norm = 0
    avg_rhf_LLL = 0
    for i in range(L.reduced_basis.nrows() - 1):
        vi = L.reduced_basis[i]
        ui = B.solve_left(vi)
        vecs_ort_q.append(ui)

    print("Gap: eta - rho = %d" % (eta - rho))

    U = Matrix(ZZ, t-1, t, vecs_ort_q)
    V = U.right_kernel()
    rec_q = V.basis()[0]
    x0 = B[0,0]

    rec_p = x0 / rec_q[0]
    print("original  p: %d" % p)
    print("recovered p: %d" % rec_p)


    if eta - 1 <= log(rec_p, 2) <= eta:
        print("ATK OK")
    else:
        print("ATK NOT OK")


def GCD_attack():

    def extract_eta_bit_factor(d):
        for fm in d.factor(): # each fm is a pair with factor and multiplicity
            f = fm[0]
            if eta-1 <= log(f, 2) <= eta:
                return f

        return False

    def test_candidate_p(candidate_p, list_samples, rho):
        for xi in list_samples:
            ri = sym_mod(xi, p)
            if abs(ri) > 2^rho:
                return False
        return True


    sec_level = 8
    rho = sec_level
    eta = 2*sec_level
    gamma = ceil((eta - rho)^2 * sec_level / (log(sec_level,2)))

    num_samples = 25

    p = sample_p(eta)

    x0 = sample_agcd(gamma, p, rho)
    x1 = sample_agcd(gamma, p, rho)

    list_samples = [sample_agcd(gamma, p, rho) for _ in range(num_samples)]

    for a in range(-2^rho, 2^rho):
        for b in range(-2^rho, 2^rho):
            d = gcd(x0 - a, x1 - b)
            if log(d, 2) >= eta-1:
                candidate_p = extract_eta_bit_factor(d)
                if candidate_p != False and test_candidate_p(candidate_p, list_samples, rho):
                    print("      d = %d" % d)
                    print("      p = %d" % p)
                    print("found p = %d" % candidate_p)
                    return True

def collision_GCD_attack():

    def remove_small_factors(a, num_fact=1000):
        q = 2
        for _ in range(num_fact):
            while q.divides(a):
                a /= q
            q = next_prime(q)
        return ZZ(a)

    sec_level = 11
    rho, eta, gamma = sec_level, 2*sec_level, 3*sec_level
    num_samples = 30
    p = sample_p(eta)

    x0 = sample_agcd(gamma, p, rho)

    mult_p = prod([x0 - r for r in range(-2^rho, 2^rho)])
    mult_p = remove_small_factors(mult_p)
    list_samples = [sample_agcd(gamma, p, rho) for _ in range(num_samples)]

    for xi in list_samples:
        tmp_mult = prod([xi - r for r in range(-2^rho, 2^rho)])
        mult_p = gcd(mult_p, tmp_mult)
        print("bitlen mult_p = %d" % mult_p.nbits())
        if eta >= log(mult_p, 2) >= eta-1:
            break


    print(" found p = %d" % mult_p)
    print(" real  p = %d" % p)
    return True


collision_GCD_attack()

   
orthogonal_lattice_atk()
