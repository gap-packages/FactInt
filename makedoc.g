##  This GAP script creates the documentation, needs: GAPDoc and AutoDoc
##  packages, pdflatex
## 
if fail = LoadPackage("AutoDoc", ">= 2016.01.21") then
    Error("AutoDoc 2016.01.21 or newer is required");
fi;

AutoDoc(rec(gapdoc:=rec(main:="main.xml")));
