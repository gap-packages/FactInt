##  This GAP script creates the documentation, needs: GAPDoc and AutoDoc
##  packages, pdflatex
## 
if fail = LoadPackage("AutoDoc", ">= 2016.01.21") then
    Error("AutoDoc 2016.01.21 or newer is required");
fi;

AutoDoc(rec( scaffold := rec(
        includes := [ "preface.xml", "general.xml", "methods.xml", "timings.xml" ],
        bib := "factintbib.xml",
        gapdoc_latex_options := rec( EarlyExtraPreamble := """
            \usepackage{amsfonts}
            \usepackage{amsxtra}
            """ ),
    )));
