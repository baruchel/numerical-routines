from sympy.matrices.dense import ones
from sympy.core import symbols
from sympy.polys.polyerrors import OptionError

def rational_interpolate(data, degnum, X=symbols('x')):
    """
    Returns a rational interpolation, where the data points are element of
    any integral domain.

    The first argument  contains the data (as a list of coordinates). The
    ``degnum`` argument is the degree in the numerator of the rational
    function. Setting it too high will decrease the maximal degree in the
    denominator for the same amount of data.

    Example:
    ========

    >>> data = [ (1,-210), (2,-35), (3,105), (4,231), (5,350), (6,465) ]
    >>> rational_interpolate(data, 2)
    (105*x**2 - 525)/(x + 1)

    Values do not need to be integers:

    >>> from sympy import sympify
    >>> x = [1, 2, 3, 4, 5, 6]
    >>> y = sympify("[-1, 0, 2, 22/5, 7, 68/7]")
    >>> rational_interpolate(zip(x,y), 2)
    (3*x**2 - 7*x + 2)/(x + 1)

    The symbol for the variable can be changed if needed:
    >>> from sympy import symbols
    >>> z = symbols('z')
    >>> rational_interpolate(data, 2, X=z)
    (105*z**2 - 525)/(z + 1)

    References
    ==========
    Algorithm is adapted from:
        http://axiom-wiki.newsynthesis.org/RationalInterpolation

    """

    xdata, ydata = list(zip(*data))

    m = degnum + 1
    k = len(xdata) - m - 1
    if k<1:
        raise OptionError("Too few values for the required degree.")
    c = ones(m+k+1, m+k+2)
    for j in range(max(m,k)):
        for i in range(m+k+1):
            c[i,j+1] = c[i,j]*xdata[i]
    for j in range(k+1):
        for i in range(m+k+1):
            c[i,m+k+1-j] = -c[i,k-j]*ydata[i]
    r = c.nullspace()[0]
    return ( sum( r[i] * X**i for i in range(m+1))
            / sum( r[i+m+1] * X**i for i in range(k+1) ) )
