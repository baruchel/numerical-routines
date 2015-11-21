# 
#   Adapted from: http://axiom-wiki.newsynthesis.org/RationalInterpolation
# 

import sympy

def rinterp(x, y, m, k):
    """Returns a rational interpolation, where the data points are element of
any integral domain.

Exactly four arguments are required: two lists containing x and f(x) values;
the degree in the expected numerator and the degree in the expected
denominator.

Code is adapted from http://axiom-wiki.newsynthesis.org/RationalInterpolation

Example:
========

    >>> from fractions import Fraction
    >>> x = [1, 2, 3, 4, 5, 6]
    >>> y=map(lambda x: sympify(Fraction(x)), ["-1","0","2","22/5","7","68/7"])
    >>> rinterp(x, y, 3, 2)
    (3*x**2 - 7*x + 2)/(x + 1)
    """
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
