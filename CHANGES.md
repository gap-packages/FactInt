This file describes changes in the FactInt package.

## 1.7.0 (unreleased)

  - Update Brent's tables of factors of `b^k +/- 1` from the collection now
    maintained by Jonathan Crombie. The full collection has grown to over
    280 MB, which is too much to distribute; the tables shipped with FactInt
    are therefore restricted to those `(b,k)` with `b <= 100`, or `b` a prime
    below 1000, or `(b,k)` already covered by FactInt 1.6.3. This is a superset
    of the data shipped with FactInt 1.6.3, adding roughly 670000 further
    factors. Use `FetchBrentFactors` to install the full collection.
  - Compress the data files in `tables/brent` in the distribution archives,
    which halves their disk usage; GAP reads them transparently

## 1.6.3 (2019-11-15)

  - Make FactInt compatible with HPC-GAP
  - Add a LICENSE file
  - Replace the internal function `PrettyInfo` by GAP's `Info` instruction,
    which improves performance in some situations
  - Minor update of Brent's tables

## 1.6.2 (2018-02-17)

  - Rewrite `FactorsTDNC` to avoid recursion
  - Optimize the loading of the data tables

## 1.6.1 (2018-01-17)

  - Update `FactorsECM` for the change that `RootInt` no longer accepts
    non-integral arguments, by converting its first argument to an integer

## 1.6.0 (2017-12-04)

  - Use Aurifeuillian factorization of `b^k + 1` for bases up to 12

## 1.5.4 (2017-02-13)

  - Rename the directory `factint/gap/` to `factint/lib/`, to follow the same
    naming convention as the subdirectories of the GAP root directory
  - Add the file `factint/doc/manual.js`

## 1.5.3 (2011-06-16)

  - Update the copy of Brent's tables of factors of integers of the form
    `b^k +/- 1`
  - Remove the CVS revision entries from the source files

## 1.5.2 (2007-09-26)

  - Use a flexible caching mechanism for factorizations of small integers,
    which yields a substantial speedup when many small numbers are factored
  - Update the copy of Brent's tables of factors of integers of the form
    `b^k +/- 1`, adding roughly 20000 new factors

## 1.5.1 (2007-09-20)

  - Convert the manual to GAPDoc format
  - Avoid unnecessarily triggering the loading of autoreadable global
    variables, which resulted in short delays and some wasting of memory

## 1.4.12 (2006-09-08)

  - Fix the formatting of the Info output of the ECM routine, which caused an
    error message if the run time spent on the first or the second stage of a
    curve was below 1ms (reported by Doug McTavish)

## 1.4.10 (2005-06-24)

  - Fix an error message when factoring a large enough integer not of the form
    `d*(10^k-1)/9` whose last 4 decimal digits were equal; present since 1.4.6
    (reported by Sven Reichard)

## 1.4.9 (2005-05-25)

  - Improve the `AbstractHTML` entry in `PackageInfo.g`

## 1.4.8 (2005-05-17)

  - Handle the special case `a^k +/- b^k`
  - Remove relics of, and compatibility with, the package loading mechanism of
    GAP 4.3

## 1.4.7 (2005-01-31)

  - Cache whole factorizations as well as single factors
  - Use the `p +/- 1` routines by default only for sufficiently large
    composites
  - Reduce the overhead for small and easy numbers further
  - Document the function `FactorsTD`
  - Add the synonyms `ECM`, `MPQS` and `CFRAC` for `FactorsECM`, `FactorsMPQS`
    and `FactorsCFRAC`

## 1.4.6 (2005-01-21)

  - Improve the performance of the factoring routine for integers of the form
    `b^k +/- 1`
  - Use Richard P. Brent's tables of factors of integers of the form
    `b^k +/- 1` (only under UNIX); the corresponding code was contributed by
    Frank Lübeck
  - Handle the following special cases:
    - two factors `p`, `q` such that `p/q` is close to a fraction with small
      numerator and denominator
    - `k! +/- 1`
    - `p1 * p2 * p3 * ... * pk +/- 1`
    - Fibonacci numbers
    - `3^k - 2^k`
    - `11111 ... 11111`
    - factors already available as values of user variables in the workspace
  - Add an option `cheap` for restricting factorization attempts to cheap
    methods

## 1.3.1 (2003-06-24)

  - Adapt the package to the new package loading and updating mechanism of
    GAP 4.4

## 1.3 (2002-07-11)

  - Reduce the overhead for factoring very small numbers
  - Add the synonym `InfoFactInt` for `IntegerFactorizationInfo`
  - Remove a restriction concerning the number of digits of some argument
    values in Info messages; this concerns mainly the ECM
  - Increase the default number of curves to be tried by `FactorsECM` for
    large numbers before invoking the MPQS
  - Introduce Info levels 4 and 5 for more frequent status messages during the
    sieving process of the MPQS

## 1.2 (2002-04-26)

  - Provide the manual also in PDF and PostScript format
  - Clean up the documentation typographically
  - Give Info messages for factorizations of small integers only if
    `InfoLevel(IntegerFactorizationInfo) = 3`, to avoid swamping users with
    messages on factorizations of numbers like 4
  - Move the declaration of `IntegerFactorizationInfo` to `factint.gd`, where
    it belongs

## 1.1

  - Use the fast special-purpose kernel routine
    `ADD_TO_LIST_ENTRIES_PLIST_RANGE` for the innermost loop of the MPQS
    sieving routine; the kernel routine was implemented in `src/listoper.c` by
    Steve Linton

## 1.0 (1999-07)

  - Initial release
