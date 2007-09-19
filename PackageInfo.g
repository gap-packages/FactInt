####################################################################################################
##
##  PackageInfo.g                         GAP4 Package `FactInt'                         Stefan Kohl
##  
#H  @(#)$Id$
##

SetPackageInfo( rec(

PackageName      := "FactInt",
Subtitle         := "Advanced Methods for Factoring Integers", 
Version          := "1.5.0",
Date             := "19/09/2007",
ArchiveURL       := "http://www.cip.mathematik.uni-stuttgart.de/~kohlsn/factint/factint-1.5.0",
ArchiveFormats   := ".tar.gz",
Persons          := [
                      rec( LastName      := "Kohl",
                           FirstNames    := "Stefan",
                           IsAuthor      := true,
                           IsMaintainer  := true,
                           Email         := "kohl@mathematik.uni-stuttgart.de",
                           WWWHome       := "http://www.cip.mathematik.uni-stuttgart.de/~kohlsn/",
                           PostalAddress := Concatenation("Institut für Geometrie und Topologie\n",
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
AbstractHTML     := Concatenation("This package provides routines for factoring integers, ",
                                  "in particular:</p>\n<ul>\n  <li>Pollard's <em>p</em>-1</li>\n",
                                  "  <li>Williams' <em>p</em>+1</li>\n  <li>Elliptic Curves ",
                                  "Method (ECM)</li>\n  <li>Continued Fraction Algorithm ",
                                  "(CFRAC)</li>\n  <li>Multiple Polynomial Quadratic Sieve ",
                                  "(MPQS)</li>\n</ul>\n<p>It also provides access to Richard P. ",
                                  "Brent's tables of factors of integers of the form ",
                                  "<em>b</em>^<em>k</em> +/- 1."),
PackageWWWHome   := "http://www.cip.mathematik.uni-stuttgart.de/~kohlsn/factint.html",
PackageDoc       := rec(
                         BookName         := "FactInt",
                         ArchiveURLSubset := ["doc"],
                         HTMLStart        := "doc/chap0.html",
                         PDFFile          := "doc/manual.pdf",
                         SixFile          := "doc/manual.six",
                         LongTitle        := "A GAP4 Package for FACToring INTegers",
                         Autoload         := true
                       ),
Dependencies     := rec(
                         GAP                    := ">=4.4.9",
                         NeededOtherPackages    := [ ["GAPDoc",">=1.0"] ],
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
                      "Continued Fraction Algorithm", "Pollard's p-1", "Williams' p+1",
                      "Cunningham Tables", "Richard P. Brent's Factor Tables" ]

) );

####################################################################################################
##
#E  PackageInfo.g  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ends here