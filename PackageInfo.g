####################################################################################################
##
##  PackageInfo.g                         GAP4 Package `FactInt'                         Stefan Kohl
##  
#H  @(#)$Id$
##

SetPackageInfo( rec(

PackageName      := "FactInt",
Subtitle         := "Advanced Methods for Factoring Integers", 
Version          := "1.4",
Date             := "04/01/2005",
ArchiveURL       := "http://www.cip.mathematik.uni-stuttgart.de/~kohlsn/factint/factint-1.4",
ArchiveFormats   := ".tar.gz",
Persons          := [
                      rec( LastName      := "Kohl",
                           FirstNames    := "Stefan",
                           IsAuthor      := true,
                           IsMaintainer  := true,
                           Email         := "kohl@mathematik.uni-stuttgart.de",
                           WWWHome       := "http://www.cip.mathematik.uni-stuttgart.de/~kohlsn/",
                           PostalAddress := Concatenation("Stefan Kohl\n",
                                                          "Institut für Geometrie und Topologie\n",
                                                          "Pfaffenwaldring 57\n",
                                                          "Universität Stuttgart\n",
                                                          "70550 Stuttgart\n",
                                                          "Germany"),
                           Place         := "Stuttgart / Germany",
                           Institution   := "University of Stuttgart"
                          )
                    ],
Status           := "accepted",
CommunicatedBy   := "Mike Atkinson (St. Andrews)",
AcceptDate       := "07/1999",
README_URL       := "http://www.cip.mathematik.uni-stuttgart.de/~kohlsn/factint/README.factint",
PackageInfoURL   := "http://www.cip.mathematik.uni-stuttgart.de/~kohlsn/factint/PackageInfo.g",
AbstractHTML     := "This package provides advanced methods for integer factorization.",
PackageWWWHome   := "http://www.cip.mathematik.uni-stuttgart.de/~kohlsn/factint.html",
PackageDoc       := rec(
                         BookName         := "FactInt",
                         ArchiveURLSubset := ["doc","htm"],
                         HTMLStart        := "htm/chapters.htm",
                         PDFFile          := "doc/manual.pdf",
                         SixFile          := "doc/manual.six",
                         LongTitle        := "A GAP4 Package for FACToring INTegers",
                         Autoload         := true
                       ),
Dependencies     := rec(
                         GAP                    := ">=4.1",
                         NeededOtherPackages    := [ ],
                         SuggestedOtherPackages := [ ],
                         ExternalConditions     := [ ]
                       ),
AvailabilityTest := ReturnTrue,
BannerString     := Concatenation( "\nLoading FactInt ", ~.Version,
                                   " (Routines for Integer Factorization )",
                                   "\nby Stefan Kohl, kohl@mathematik.uni-stuttgart.de\n\n" ),
Autoload         := true,
TestFile         := "factint.tst",
Keywords         := [ "Integer factorization", "ECM", "Elliptic Curves Method",
                      "MPQS", "Multiple Polynomial Quadratic Sieve", "CFRAC",
                      "Continued Fraction Algorithm", "Pollard's p-1", "Williams' p+1" ]

) );

####################################################################################################
##
#E  PackageInfo.g  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ends here