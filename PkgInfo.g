########################################################################################################################################################################
##
##  PkgInfo.g                                                             GAP4 Package `FactInt'                                                             Stefan Kohl
##  
#H  @(#)$Id$
##

SetPackageInfo( rec(

PkgName          := "FactInt",
Version          := "1.2",
Date             := "30/04/2002",
ArchiveURL       := "http://www.cip.mathematik.uni-stuttgart.de/~kohlsn/factint/factint-1.2",
ArchiveFormats   := ".zoo",
Persons          := [
                      rec( LastName      := "Kohl",
                           FirstNames    := "Stefan",
                           IsAuthor      := true,
                           IsMaintainer  := true,
                           Email         := "kohl@mathematik.uni-stuttgart.de",
                           WWWHome       := "http://www.cip.mathematik.uni-stuttgart.de/~kohlsn/",
                           PostalAddress := "Stefan Kohl\nMathematisches Institut B, 2. Lehrstuhl\nPfaffenwaldring 57\nUniversität Stuttgart\n70550 Stuttgart\nGermany",
                           Place         := "Stuttgart / Germany",
                           Institution   := "University of Stuttgart"
                          )
                    ],
Status           := "accepted",
CommunicatedBy   := "Mike Atkinson (St. Andrews)",
AcceptDate       := "07/1999",
README_URL       := "http://www.cip.mathematik.uni-stuttgart.de/~kohlsn/factint/README.factint",
PkgInfoURL       := "http://www.cip.mathematik.uni-stuttgart.de/~kohlsn/factint/PkgInfo.g",
AbstractHTML     := "This package provides advanced methods for integer factorization.",
PackageWWWHome   := "http://www.cip.mathematik.uni-stuttgart.de/~kohlsn/factint.html",
PackageDoc       := rec(
                         BookName  := "FactInt",
                         Archive   := "http://www.cip.mathematik.uni-stuttgart.de/~kohlsn/factint/factint-1.2doc-win.zip",
                         HTMLStart := "htm/chapters.htm",
                         PDFFile   := "doc/manual.pdf",
                         SixFile   := "doc/manual.six",
                         LongTitle := "A GAP4 Package for FACToring INTegers",
                         AutoLoad  := true
                       ),
Dependencies     := rec(
                         GAP                    := ">=4.1",
                         NeededOtherPackages    := [ ],
                         SuggestedOtherPackages := [ ],
                         ExternalConditions     := [ ]
                       ),
AvailabilityTest := ReturnTrue,
Autoload         := true,
TestFile         := "factint.tst",
Keywords         := [ "Integer factorization", "ECM", "Elliptic Curves Method",
                      "MPQS", "Multiple Polynomial Quadratic Sieve", "CFRAC",
                      "Continued Fraction Algorithm", "Pollard's p-1", "Williams' p+1" ]

) );

########################################################################################################################################################################
##
#E  PkgInfo.g  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ends here
