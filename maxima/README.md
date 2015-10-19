# Numerical routines for Maxima

### rinterp

Rational interpolation (adapted from http://axiom-wiki.newsynthesis.org/RationalInterpolation ).

    (%i1) xlist: [0, 1/2, 3/4, 2, 7/3]$
    (%i2) ylist: map( lambda([x], (2*x^2-5*x+7)/(4*x^2-5)), xlist)$
    (%i3) rinterp( xlist, ylist, 2, 2);
                                      2
                                     x    5 x   7
                                     -- - --- + -
                                     2     4    4
    (%o3)                            ------------
                                         2   5
                                        x  - -
                                             4
