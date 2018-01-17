####################################################################################################
##
##  PackageInfo.g                         GAP4 Package `FactInt'                         Stefan Kohl
##  
####################################################################################################

SetPackageInfo( rec(

PackageName      := "FactInt",
Subtitle         := "Advanced Methods for Factoring Integers", 
Version          := "1.6.0",
Date             := "04/12/2017",
ArchiveURL       := "https://stefan-kohl.github.io/factint/factint-1.6.0",
ArchiveFormats   := ".tar.gz", # "-win.zip" when providing text files with Windows line breaks
Persons          := [
                      rec( LastName      := "Kohl",
                           FirstNames    := "Stefan",
                           IsAuthor      := true,
                           IsMaintainer  := true,
                           Email         := "stefan@gap-system.org",
                           WWWHome       := "https://stefan-kohl.github.io/"
                         )
                    ],
Status           := "accepted",
CommunicatedBy   := "Mike Atkinson (St. Andrews)",
AcceptDate       := "07/1999",
PackageWWWHome   := "https://stefan-kohl.github.io/factint.html",
README_URL       := "https://stefan-kohl.github.io/factint/README.factint",
PackageInfoURL   := "https://stefan-kohl.github.io/factint/PackageInfo.g",
AbstractHTML     := Concatenation("This package provides routines for factoring integers, ",
                                  "in particular:</p>\n<ul>\n  <li>Pollard's <em>p</em>-1</li>\n",
                                  "  <li>Williams' <em>p</em>+1</li>\n  <li>Elliptic Curves ",
                                  "Method (ECM)</li>\n  <li>Continued Fraction Algorithm ",
                                  "(CFRAC)</li>\n  <li>Multiple Polynomial Quadratic Sieve ",
                                  "(MPQS)</li>\n</ul>\n<p>It also provides access to Richard P. ",
                                  "Brent's tables of factors of integers of the form ",
                                  "<em>b</em>^<em>k</em> +/- 1."),
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
                         GAP                    := ">=4.8.8",
                         NeededOtherPackages    := [ ["GAPDoc",">=1.6"] ],
                         SuggestedOtherPackages := [ ],
                         ExternalConditions     := [ ]
                       ),
AvailabilityTest := ReturnTrue,
BannerString     := Concatenation( "\nLoading FactInt ", ~.Version,
                                   " (Routines for Integer Factorization)",
                                   "\nby Stefan Kohl, stefan@gap-system.org\n\n" ),
TestFile         := "tst/testall.g",
Keywords         := [ "Integer factorization", "ECM", "Elliptic Curves Method",
                      "MPQS", "Multiple Polynomial Quadratic Sieve", "CFRAC",
                      "Continued Fraction Algorithm", "Pollard's p-1", "Williams' p+1",
                      "Cunningham Tables", "Richard P. Brent's Factor Tables" ],

AutoDoc := rec(
  TitlePage := rec(
    Copyright := """
      &copyright; 1999 - 2017 by Stefan Kohl. <P/>

      &FactInt; is free software: you can redistribute it and/or modify
      it under the terms of the GNU General Public License as published by
      the Free Software Foundation, either version 2 of the License, or
      (at your option) any later version. <P/>

      &FactInt; is distributed in the hope that it will be useful,
      but WITHOUT ANY WARRANTY; without even the implied warranty of
      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
      GNU General Public License for more details. <P/>

      For a copy of the GNU General Public License, see
      the file <F>GPL</F> in the <F>etc</F> directory of the &GAP;
      distribution or see <URL>http://www.gnu.org/licenses/gpl.html</URL>.
      """,
    Abstract := """<#Include SYSTEM "abstract.xml">""",
    Acknowledgements := """
      I would like to thank Bettina Eick and Steve Linton for their support
      and many interesting discussions.
      """,
  ),
),

) );

####################################################################################################
##
#E  PackageInfo.g  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ends here
