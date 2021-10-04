Zx.<x> = ZZ['x']

def inv_g_ZZ(z, B, l):
    return vector(ZZ(z).digits(base=B, padto=l))

def inv_g_poly(poly, B, l, n):
    result = vector(Zx, [0] * l)
    poly = Zx(poly.list())
    pow_x = 1
    for i in range(n):
        result += pow_x * inv_g_ZZ(poly[i], B, l)
        pow_x *= x
    return result
