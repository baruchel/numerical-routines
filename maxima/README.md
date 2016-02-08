# Numerical routines for Maxima

### rinterp.mac

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

### pslq.mac

An implementation of the PSLQ algorithm for detecting integer relations between constants.

    (%i1) load("pslq.mac");
    (%o1)                              pslq.mac
    (%i2) fpprec: 72$
    (%i3) a:bfloat(5/4*%pi-7/3*%e+1/4)$
    (%i4) pslq([a, %pi, %e, 1],1b-64);
    (%o4)                         [12, - 15, 28, - 3]
    (%i5) float(12*a - 15*%pi + 28*%e - 3);
    (%o5)                                 0.0

### convolution.lisp

Low-level Lisp code for various operations related to convolutions of sequences. A very fast function for computing the generating function is provided (see the content of the file for explanations on how to compile the file and use it).
