README file for `FactInt' package.

This package for GAP 4 provides routines for integer factorization,
in particular:

 - Pollard's p-1
 - Williams' p+1
 - The Elliptic Curves Method (ECM)
 - The Continued Fraction Algorithm (CFRAC)
 - The Multiple Polynomial Quadratic Sieve (MPQS)

It also provides access to Richard P. Brent's tables of factors of
integers of the form b^k +/- 1.

FactInt is completely written in the GAP language. It neither contains
nor requires external binaries. No other packages are needed.
It must be installed in the pkg subdirectory of the GAP distribution.

By default, `FactInt' is autoloaded.
You can load the package also via LoadPackage( "factint" );

If you have problems with this package, wish to make comments
or suggestions, or if you find bugs, please send e-mail to

kohl@mathematik.uni-stuttgart.de
