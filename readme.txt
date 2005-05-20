

                            The FactInt Package
                            ===================


                                 Abstract

FactInt is a GAP 4 package which provides routines for factoring integers,
in particular:

 - Pollard's p-1
 - Williams' p+1
 - The Elliptic Curves Method (ECM)
 - The Continued Fraction Algorithm (CFRAC)
 - The Multiple Polynomial Quadratic Sieve (MPQS)

It also provides access to Richard P. Brent's tables of factors of integers
of the form b^k +/- 1.


                               Requirements

The FactInt Package needs at least version 4.4 of GAP, is completely written
in the GAP language and does neither contain nor require external binaries.
No other packages are needed.


                               Installation

Like any other GAP package, FactInt must be installed in the pkg/
subdirectory of the GAP distribution. This is accomplished by extracting
the distribution file in this directory. By default, FactInt is
autoloaded. This means that it is loaded automatically when you start GAP.


If you have problems with this package, wish to make comments or
suggestions, or if you find bugs, please send e-mail to

Stefan Kohl, kohl@mathematik.uni-stuttgart.de

