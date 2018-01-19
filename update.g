#
# GitHubPagesForGAP - a template for using GitHub Pages within GAP packages
#
# Copyright (c) 2013-2014 Max Horn
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#

# Parse PackageInfo.g and regenerate _data/package.yml from it.

PrintPeopleList := function(stream, people)
    local p;
    for p in people do
        AppendTo(stream, "    - name: ", p.FirstNames, " ", p.LastName, "\n");
        if IsBound(p.WWWHome) then
            AppendTo(stream, "      url: ", p.WWWHome, "\n");
        elif IsBound(p.Email) then
            AppendTo(stream, "      url: mailto:", p.Email, "\n");
        fi;
    od;
    AppendTo(stream, "\n");
end;

PrintPackageList := function(stream, pkgs)
    local p, pkginfo;
    for p in pkgs do
        AppendTo(stream, "    - name: \"", p[1], "\"\n");
        AppendTo(stream, "      version: \"", p[2], "\"\n");
        pkginfo := PackageInfo(p[1]);
        if Length(pkginfo) > 0 and IsBound(pkginfo[1].PackageWWWHome) then
            AppendTo(stream, "      url: \"", pkginfo[1].PackageWWWHome, "\"\n");
        fi;
    od;
    AppendTo(stream, "\n");
end;

GeneratePackageYML:=function(pkg)
    local stream, authors, maintainers, contributors, formats, f, tmp;

    stream := OutputTextFile("_data/package.yml", false);
    SetPrintFormattingStatus(stream, false);
    
    AppendTo(stream, "name: ", pkg.PackageName, "\n");
    AppendTo(stream, "version: ", pkg.Version, "\n");
    AppendTo(stream, "date: ", pkg.Date, "\n"); # TODO: convert to ISO 8601?
    AppendTo(stream, "description: |\n");
    AppendTo(stream, "    ", pkg.Subtitle, "\n");
    AppendTo(stream, "\n");

    authors := Filtered(pkg.Persons, p -> p.IsAuthor);
    if Length(authors) > 0 then
        AppendTo(stream, "authors:\n");
        PrintPeopleList(stream, authors);
    fi;

    maintainers := Filtered(pkg.Persons, p -> p.IsMaintainer);
    if Length(maintainers) > 0 then
        AppendTo(stream, "maintainers:\n");
        PrintPeopleList(stream, maintainers);
    fi;

    contributors := Filtered(pkg.Persons, p -> not p.IsMaintainer and not p.IsAuthor);
    if Length(contributors) > 0 then
        AppendTo(stream, "contributors:\n");
        PrintPeopleList(stream, contributors);
    fi;

    if IsBound(pkg.Dependencies.GAP) then
        AppendTo(stream, "GAP: \"", pkg.Dependencies.GAP, "\"\n\n");
    fi;

    if IsBound(pkg.Dependencies.NeededOtherPackages) and
        Length(pkg.Dependencies.NeededOtherPackages) > 0 then
        AppendTo(stream, "needed-pkgs:\n");
        PrintPackageList(stream, pkg.Dependencies.NeededOtherPackages);
    fi;

    if IsBound(pkg.Dependencies.SuggestedOtherPackages) and
        Length(pkg.Dependencies.SuggestedOtherPackages) > 0 then
        AppendTo(stream, "suggested-pkgs:\n");
        PrintPackageList(stream, pkg.Dependencies.SuggestedOtherPackages);
    fi;

    AppendTo(stream, "www: ", pkg.PackageWWWHome, "\n");
    tmp := SplitString(pkg.README_URL,"/");
    tmp := tmp[Length(tmp)];  # extract README filename (typically "README" or "README.md")
    AppendTo(stream, "readme: ", tmp, "\n");
    AppendTo(stream, "packageinfo: ", pkg.PackageInfoURL, "\n");
    if IsBound(pkg.GithubWWW) then
        AppendTo(stream, "github: ", pkg.GithubWWW, "\n");
    fi;
    AppendTo(stream, "\n");
    
    formats := SplitString(pkg.ArchiveFormats, " ");
    if Length(formats) > 0 then
        AppendTo(stream, "downloads:\n");
        for f in formats do
            AppendTo(stream, "    - name: ", f, "\n");
            AppendTo(stream, "      url: ", pkg.ArchiveURL, f, "\n");
        od;
        AppendTo(stream, "\n");
    fi;

    AppendTo(stream, "abstract: |\n");
    for tmp in SplitString(pkg.AbstractHTML,"\n") do
        AppendTo(stream, "    ", tmp, "\n");
    od;
    AppendTo(stream, "\n");

    AppendTo(stream, "status: ", pkg.Status, "\n");
    if IsRecord(pkg.PackageDoc) then
        AppendTo(stream, "doc-html: ", pkg.PackageDoc.HTMLStart, "\n");
        AppendTo(stream, "doc-pdf: ", pkg.PackageDoc.PDFFile, "\n");
    else
        Assert(0, IsList(pkg.PackageDoc));
        AppendTo(stream, "doc-html: ", pkg.PackageDoc[1].HTMLStart, "\n");
        AppendTo(stream, "doc-pdf: ", pkg.PackageDoc[1].PDFFile, "\n");
        if Length(pkg.PackageDoc) > 1 then
            Print("Warning, this package has more than one help book!\n");
        fi;
    fi;

    # TODO: use Keywords?

    CloseStream(stream);
end;
Read("PackageInfo.g");
GeneratePackageYML(GAPInfo.PackageInfoCurrent);
QUIT;
