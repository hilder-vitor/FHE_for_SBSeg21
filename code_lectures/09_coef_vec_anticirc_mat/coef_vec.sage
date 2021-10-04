load("utils.sage")

#   Given a polynomial g, this function returns a vector representing g % f.
#   The dimension of the returned polynomial is deg(f).
#   Notice that this map is an isomorphism from ZZ[x]/<x^N+1> to ZZ^N.
def poly_to_vec(g, N):
    f = x^N+1
    v = (g % f).coefficients(sparse=False)
    return vector(v + [0]*(N - len(v)))

#   Given a polynomial g, this function returns a matrix that represents
# the product by g in the ring ZZ[x]/<x^N+1>, that is,
# a matrix A that satisfies A * poly_to_vec(h) = poly_to_vec(h*g % f) for all h.
#   The returned matrix is NxN.
def poly_to_mat(g, N):
    A = Matrix(ZZ, [poly_to_vec(g * x^i, N) for i in range(N)])
    #return A.transpose()
    return A
