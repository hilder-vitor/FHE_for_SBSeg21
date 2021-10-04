Zx.<x> = ZZ['x']

def sym_mod(a, n):
    a = ZZ(a) % n
    if 2*a > n:
        return a - n
    return a

def sym_mod_poly(poly, q):
    return Zx([sym_mod(ZZ(ai), q) for ai in poly.list()])

def sym_mod_vec(vec, q):
    return [sym_mod_poly(vi, q) for vi in vec]

def round_poly(h):
    return Zx([round(hi) for hi in h.list()])

def infinity_norm(poly):
    if 0 == poly:
        return 0
    else:
        return vector(ZZ, poly.list()).norm(Infinity)

def infinity_norm_vec(vec):
    return max([infinity_norm(vi) for vi in vec])

def random_bit():
    if random() > 0.5:
        return 1
    return 0


