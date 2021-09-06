def sample_r(rho):
    return ZZ.random_element(-2^rho,2^rho)
def sample_q(gamma, eta):
    return ZZ.random_element(0, 2^(gamma - eta))
def sample_acd(gamma, p, rho):
    eta = ceil(log(p, 2))
    r = sample_r(rho)
    q = sample_q(gamma, eta)
    return p*q + r
