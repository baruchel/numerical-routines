# Numerical routines for Pari-GP

### rinterp

Rational interpolation (adapted from http://axiom-wiki.newsynthesis.org/RationalInterpolation ).

    ? xlist = [0, 1/2, 3/4, 2, 7/3]
    %1 = [0, 1/2, 3/4, 2, 7/3]
    ? ylist = vector(length(xlist), k, x=xlist[k]; (2*x^2-5*x+7)/(4*x^2-5))
    %2 = [-7/5, -5/4, -35/22, 5/11, 56/151]
    ? rinterp(xlist, ylist, 2, 2)
    %3 = (2*x^2 - 5*x + 7)/(4*x^2 - 5)
