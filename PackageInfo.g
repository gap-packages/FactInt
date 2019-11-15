####################################################################################################
##
##  PackageInfo.g                         GAP4 Package `FactInt'                         Stefan Kohl
##  
####################################################################################################

SetPackageInfo( rec(

PackageName      := "FactInt",
Subtitle         := "Advanced Methods for Factoring Integers", 
Version          := "1.6.3",
Date             := "15/11/2019", # dd/mm/yyyy format
License          := "GPL-2.0-or-later",

SourceRepository := rec(
    Type := "git",
    URL := Concatenation( "https://github.com/gap-packages/", ~.PackageName ),
),
IssueTrackerURL := Concatenation( ~.SourceRepository.URL, "/issues" ),
PackageWWWHome  := Concatenation( "https://gap-packages.github.io/", ~.PackageName ),
README_URL      := Concatenation( ~.PackageWWWHome, "/README.md" ),
PackageInfoURL  := Concatenation( ~.PackageWWWHome, "/PackageInfo.g" ),
ArchiveURL      := Concatenation( ~.SourceRepository.URL,
                                 "/releases/download/v", ~.Version,
                                 "/", ~.PackageName, "-", ~.Version ),
ArchiveFormats := ".tar.gz",

Persons          := [
                      rec( LastName      := "Kohl",
                           FirstNames    := "Stefan",
                           IsAuthor      := true,
                           IsMaintainer  := true,
                           Email         := "stefan@gap-system.org",
                           WWWHome       := "https://stefan-kohl.github.io/"
                         ),
                      rec( LastName      := "Konovalov",
                           FirstNames    := "Alexander",
                           IsAuthor      := false,
                           IsMaintainer  := true,
                           Email         := "alexander.konovalov@st-andrews.ac.uk",
                           WWWHome       := "https://alexk.host.cs.st-andrews.ac.uk",
                           PostalAddress := Concatenation( [
                             "School of Computer Science\n",
                             "University of St Andrews\n",
                             "Jack Cole Building, North Haugh,\n",
                             "St Andrews, Fife, KY16 9SX, Scotland" ] ),
                           Place         := "St Andrews",
                           Institution   := "University of St Andrews"
    ),
                    ],
Status           := "accepted",
CommunicatedBy   := "Mike Atkinson (St. Andrews)",
AcceptDate       := "07/1999",
AbstractHTML     := """
This package provides routines for factoring integers, in particular:
<ul>
  <li>Pollard's <em>p</em>-1</li>
  <li>Williams' <em>p</em>+1</li>
  <li>Elliptic Curves Method (ECM)</li>
  <li>Continued Fraction Algorithm (CFRAC)</li>
  <li>Multiple Polynomial Quadratic Sieve (MPQS)</li>
</ul>
It also provides access to Richard P. Brent's tables of factors of integers of the form <em>b</em>^<em>k</em> +/- 1.
""",

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
