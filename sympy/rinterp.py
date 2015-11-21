# 
#   Adapted from: http://axiom-wiki.newsynthesis.org/RationalInterpolation
# 

import sympy

def rinterp(x, y, m, k):
    if len(x) != len(y):
        raise sympy.ShapeError("Different number of points and values.")
    if len(x) != m+k+1:
        raise sympy.ShapeError("Wrong number of points.")
    c = sympy.ones(m+k+1, m+k+2)
    for j in range(max(m,k)):
        for i in range(m+k+1):
            c[i,j+1] = c[i,j]*x[i]
    for j in range(k+1):
        for i in range(m+k+1):
            c[i,m+k+1-j] = -c[i,k+1-j]*y[i]
    r = c.nullspace()[0]
    z = sympy.var('x')
    return ( sum( r[i+1] * z**i for i in range(m))
            / sum( r[i+m+1] * z**i for i in range(k) ) )
