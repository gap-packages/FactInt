#############################################################################
##
#W  general.gi              GAP4 Package `FactInt'                Stefan Kohl
##                                                               Frank Lübeck
##
#H  @(#)$Id$
##
##  This file contains the general routines for integer factorization and
##  auxiliary functions used by them and/or more than one of the 
##  functions for the specific factorization methods implemented in
##  pminus1.gi (Pollard's $p-1$), pplus1.gi (Williams' $p+1$), ecm.gi 
##  (Elliptic Curves Method, ECM), cfrac.gi (Continued Fraction Algorithm,
##  CFRAC) and mpqs.gi (Multiple Polynomial Quadratic Sieve, MPQS).
##
##  In each algorithm, <n> is the number to be factored.
## 
##  Descriptions of the algorithms can be found in
##
##  David M. Bressoud: Factorization and Primality Testing, Springer 1989
##
##  A (brief) description of the factoring algorithms can also be found in
##
##  Henri Cohen: A Course in Computational Algebraic Number Theory,
##  Springer 1993
##
##  In the last book, there is also a (very short) description of the
##  Generalized Number Field Sieve (GNFS), which is the most efficient
##  factoring method known today, but is not implemented here, because
##  the MPQS is usually faster for numbers less than $10^{100}$, say,
##  and factoring ``difficult'' numbers of this order of magnitude 
##  is far beyond the scope in this context.
##
Revision.general_gi :=
  "@(#)$Id$";

InstallGlobalFunction( FactIntInfo,
                       function( lev ) 
                         SetInfoLevel(IntegerFactorizationInfo,lev); 
                       end );


# For pretty-printing info messages

PrettyInfo := function (lev,Args)

  local InfoString,Arg;
  
  InfoString := "";
  for Arg in Args do
    if   IsString(Arg) 
    then Append(InfoString,Arg);
    else Append(InfoString,String(Arg[1],Arg[2]));
    fi;
  od;
  Info(IntegerFactorizationInfo,lev,InfoString);
end;
MakeReadOnlyGlobal("PrettyInfo");


# For converting a time in ms as given by Runtime() to a
# printable string

TimeToString := function (Time)

  return Concatenation(String(Int(Time/1000)),".",
                       String(Time mod 1000 + 1000){[2..4]}," sec.");
end;
MakeReadOnlyGlobal("TimeToString");


# For checking the results of all the factorization routines

FactorizationCheck := function (n,Result)

  local  ResultCorrect;

  if IsList(Result[1])
  then ResultCorrect :=     Product(Flat(Result)) = n 
                        and ForAll(Result[1],IsProbablyPrimeInt)
                        and not ForAny(Result[2],IsProbablyPrimeInt);
  else ResultCorrect :=     Product(Result) = n
                        and ForAll(Result,IsProbablyPrimeInt);
  fi;
  if not ResultCorrect
  then Error("\nInternal error, the result is incorrect !!!\n\n",
             "Please send e-mail to the author\n",
             "(kohl@mathematik.uni-stuttgart.de)\n",
             "and mention the number to be factored : \n",n,
             "\nas well as the options you specified, ",
             "thank you very much.\n"); 
  fi;
end;
MakeReadOnlyGlobal("FactorizationCheck");


# For writing the temporary factorization data of the MPQS
# (relations over the factor base etc.) to a file which can
# be read using the `Read'-function

SaveMPQSTmp := function (TempFile)

  local  MPQSTmp;

  MPQSTmp := ValueOption("MPQSTmp");
  PrintTo(TempFile,
          "PushOptions(rec(MPQSTmp :=\n",MPQSTmp,"));\n");
end;
MakeReadOnlyGlobal("SaveMPQSTmp");


# Initialize the prime differences list
# (used by ECM, Pollard's $p-1$ and Williams' $p+1$ for second stages)

BindGlobal("PrimeDiffs",[]);
BindGlobal("PrimeDiffLimit",1000000);

InitPrimeDiffs := function (Limit)

  local  Sieve,p,Maxp,pos,incr,zero,one;

  if Limit <= PrimeDiffLimit and PrimeDiffs <> [] 
  then return; fi;
  Limit := Maximum(Limit,PrimeDiffLimit);
  Info(IntegerFactorizationInfo,2,
       "Initializing prime differences list, ",
       "PrimeDiffLimit = ",Limit);
  MakeReadWriteGlobal("PrimeDiffLimit");
  PrimeDiffLimit := Limit;
  MakeReadOnlyGlobal("PrimeDiffLimit"); 
  zero := Zero(GF(2)); one := One(GF(2));
  Sieve := ListWithIdenticalEntries(PrimeDiffLimit,zero);
  Sieve[1] := one;
  Maxp := RootInt(PrimeDiffLimit); p := 2;
  while p <= Maxp do
    pos := 2 * p;
    while pos <= PrimeDiffLimit do
      Sieve[pos] := one;
      pos := pos + p;
    od;
    p := NextPrimeInt(p);  
  od;
  MakeReadWriteGlobal("PrimeDiffs");
  PrimeDiffs := [2,1];
  incr := 0;
  for pos in [4..PrimeDiffLimit] do
    incr := incr + 1;
    if Sieve[pos] = zero then
      Add(PrimeDiffs,incr);
      incr := 0;
    fi;
  od;
  MakeReadOnlyGlobal("PrimeDiffs");
end;
MakeReadOnlyGlobal("InitPrimeDiffs");


# BRENTFACTORS is a list of lists. If there is an entry in position [a][n]
# then this is a list of primes which divide b^k - 1 but no b^l - 1 with 
# l < k.
#  
# BRENTFACTORSAVAILABLE is a binary list (blist) whose entry [b] is 'true'
# when there are data available for BRENTFACTORS[b]. (They are only loaded 
# when they are needed.)
#
# The source for the data is
#
# http://web.comlab.ox.ac.uk/oucl/work/richard.brent/factors.html.
#
# The code for accessing these factorization tables has been contributed
# by Frank Lübeck.

BindGlobal("BRENTFACTORS", []);
BindGlobal("BRENTFACTORSAVAILABLE", BlistStringDecode(
"6E7EFF5EEFFF7FFEFFFF7FFFEFFFFF76FFFEFFFFFFFFFFFFEFFFFEFFFFFFDFFEFFFFFFFFFFFFF\
FFFFFFFFDFFFFFFFFFFFFFFFFFFFFFFFF7FFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFF\
F7FFFFFFFFFFFFFFFFFFFFFFFFF7FFFFFFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFF\
FFFFFFF7FFFFFFFFE0084042408s0202s02032800101A60s02020020s03080881s03010D820022\
200014s0210000100A0200008s032004s033024s0202s04280240s0301s05080402s0280s03088\
0s0220s0302s020210s080822s0701s0210s0780s0980s0320s02A0s0708000210s02200002s02\
02s03082000080002s0A20s07040082s078080s0B1040s07020080s0520s2208s02040080s0518\
s04200240s0740s0980s0C02s0480s0320s0702s0C408002s0402s0320s0740s09040040s08800\
020s0240s0C02s0308s1680s110201s0620s0320s050840s1708s0A02s0B82s0D0280s0B02s0A2\
0s0A20s0808s1008s0D20s1E80s0F20s2880s1E08s5880s1C02s0302s0D2008s4602s0902s3208\
s1508s0208s0508s1120"
));

BindGlobal( "WriteBrentFactorsFiles",

  function ( dir )

    local  bf, i;

    if not IsDirectory(dir) then dir := Directory(dir); fi;
    bf := BRENTFACTORS;
    for i in [1..Length(bf)] do 
      if IsBound(bf[i]) then 
        PrintTo(Filename(dir, Concatenation("brfac", String(i))),
                "BRENTFACTORS[",i,"]:=", bf[i], ";\n"); 
      fi;
    od;
  end );

#############################################################################
##
#F  FetchBrentFactors( ) . . get Brent's tables of factors of numbers b^k - 1
##
InstallGlobalFunction( "FetchBrentFactors",

  function ( )

    local  str, get, comm, rows, b, k, a, dir;

    # fetch the file from R. Brents ftp site and gunzip it into 'str'

    str := "";
    get := OutputTextString(str, false);
    comm := Concatenation("wget -q ftp://ftp.comlab.ox.ac.uk/pub/",
                          "Documents/techpapers/Richard.Brent/factors/",
                          "factors.gz -O - | gzip -dc ");
    Process(DirectoryCurrent(), Filename(DirectoriesSystemPrograms(),"sh"),
            InputTextUser(), get, ["-c", comm]);
  
    rows := SplitString(str, "", "\n");
    str := 0;
    for a in rows do 
      b := List(SplitString(a, "", "+- \n"), Int);
      if not IsBound(BRENTFACTORS[b[1]]) then
        BRENTFACTORS[b[1]] := [];
      fi;
      if '-' in a then
        k := b[2];
      else
        k := 2*b[2];
      fi;
      if not IsBound(BRENTFACTORS[b[1]][k]) then
        BRENTFACTORS[b[1]][k] := [b[3]];
      else
        Add(BRENTFACTORS[b[1]][k], b[3]);
      fi;
    od;
    dir := GAPInfo.PackagesInfo.("factint")[1].InstallationPath;
    WriteBrentFactorsFiles(Concatenation(dir,"/tables/"));
  end );


# Grab factors with at least <mindigits> decimal digits from <file>.
# Optionally exclude rightmost (largest?) factor.

GrabFactors := function ( file, mindigits, excludelast )

  local  nums, num, lines, line, nondigits;

  lines     := SplitString(ReadAll(InputTextFile(file)),"\n","\n");
  nondigits := Difference(List([0..255],CHAR_INT),"0123456789");
  nums  := [];
  for line in lines do
    num := List(SplitString(line,nondigits,nondigits),Int);
    if IsEmpty(num) then continue; fi;
    if excludelast then Unbind(num[Length(num)]); fi;
    num := Filtered(num,n->LogInt(n,10)>=mindigits-1);
    num := Filtered(num,IsProbablyPrimeInt);
    nums := Concatenation(nums,num);
  od;
  return Set(nums);  
end;
MakeReadOnlyGlobal("GrabFactors");

# Remove redundancies from a list <facts> of factors of numbers <f>(k)
# for 1 <= k <= <max_k>, under the assumption that a|b implies f(a)|f(b).

CleanedFactorsList := function ( facts, f, max_k )

  local  result, smldivpos, val, i;

  val := List([1..max_k],f);
  smldivpos := List(facts,p->First([1..max_k],k->val[k] mod p = 0));
  result := List([1..max_k],k->facts{Filtered([1..Length(smldivpos)],
                                              i->smldivpos[i]=k)});;
  for i in [1..Length(result)] do
    if result[i] <> [] then Unbind(result[i][Length(result[i])]); fi;
  od;
  return Set(Flat(result));
end;
MakeReadOnlyGlobal("CleanedFactorsList");


# Apply a factoring method to the composite factors of a partial
# factorization and give information about it

ApplyFactoringMethod := function (arg)

  local  FactoringMethod,Parameters,FactList,Bound,
         InfoArgs,InfoArgsTmp,InfoBaseString,InfoString,Display_n,
         Unfactored,n,Arguments,Temp,l;

  FactoringMethod := arg[1];
  Parameters      := arg[2];
  FactList        := arg[3];
  Bound           := arg[4];
  if FactList[2] = [] then return; fi;
  Display_n := false;
  if IsBound(arg[5]) then
    InfoArgs := arg[5];
    l := Length(InfoArgs);
    Display_n := InfoArgs[l] = "n";
    if Display_n then Unbind(InfoArgs[l]); fi;
    InfoArgsTmp := List(InfoArgs,function(Arg) 
                                   if IsFunction(Arg) or Arg = fail
                                   then return "<func.>";
                                   else return Arg; fi;
                                 end);
    InfoBaseString := Concatenation(List(InfoArgsTmp,elt->String(elt))); 
    if not Display_n 
    then Info(IntegerFactorizationInfo,2,"");
         Info(IntegerFactorizationInfo,1,InfoBaseString); fi;
  fi;
  Unfactored := ShallowCopy(FactList[2]);
  FactList[2] := [];
  for n in Unfactored do
    if n < Bound then
      if Display_n then 
        if   Length(InfoBaseString) 
           + LogInt(n,10) + 1 >= SizeScreen()[1]
        then InfoString := Concatenation(InfoBaseString,"\n");
        else InfoString := InfoBaseString; fi;
        Info(IntegerFactorizationInfo,2,"");
        Info(IntegerFactorizationInfo,1,
             Concatenation(InfoString,String(n)));
      fi;
      Arguments := [n];
      Append(Arguments,Parameters);
      Temp := CallFuncList(FactoringMethod,Arguments);
      if not IsList(Temp[1]) then Temp := [Temp,[]]; fi;
      if Temp[1] <> [] then
      Info(IntegerFactorizationInfo,1,"Intermediate result : ",Temp); fi;
      Append(FactList[1],Temp[1]); Sort(FactList[1]);
      Append(FactList[2],Temp[2]); Sort(FactList[2]);
    else Add(FactList[2],n);
    fi;
  od;
end;
MakeReadOnlyGlobal("ApplyFactoringMethod");


# Trial Division

FactorsTD := function (arg)

  local n,p,Result,DivisorsList;

  n := arg[1];
  if IsBound(arg[2]) then DivisorsList := arg[2];
                     else DivisorsList := Primes; fi;
  Result := [[],[]];
  for p in DivisorsList do
    while n mod p = 0 do 
      if IsProbablyPrimeInt(p) then Add(Result[1],p); 
                       else Add(Result[2],p); fi;
      n := n/p;
      if IsProbablyPrimeInt(n) then Add(Result[1],n); n := 1; fi;
    od;
    if n = 1 then return Result; fi;
  od;
  if IsProbablyPrimeInt(n) then Add(Result[1],n);
                           else Add(Result[2],n); fi;   
  return Result;
end;
MakeReadOnlyGlobal("FactorsTD");

# Initialize some lists of trial divisors

BindGlobal("K_FACTORIAL_M1_FACTORS",[]);
BindGlobal("K_FACTORIAL_P1_FACTORS",[]);
BindGlobal("K_PRIMORIAL_M1_FACTORS",[]);
BindGlobal("K_PRIMORIAL_P1_FACTORS",[]);
BindGlobal("FACTORS_FIB",[]);
BindGlobal("FIB_RES", # Fib(k) mod 13, 21, 34, 55, 89, 144.
[ [ 0, 1, 2, 3, 5, 8, 10, 11, 12 ], [ 0, 1, 2, 3, 5, 8, 13, 18, 20 ],
  [ 0, 1, 2, 3, 5, 8, 13, 21, 26, 29, 31, 32, 33 ],
  [ 0, 1, 2, 3, 5, 8, 13, 21, 34, 47, 52, 54 ],
  [ 0, 1, 2, 3, 5, 8, 13, 21, 34, 55, 68, 76, 81, 84, 86, 87, 88 ],
  [ 0, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 123, 136, 141, 143 ] ]);
BindGlobal("POW3_M_POW2_FACTORS",[]);


# Treat values of functions f such that a|b implies f(a)|f(b)
# (f is assumed to be strictly growing)

FactorsMultFunc := function ( n, f )

  local  val, fact, k, step, fk, gcd, i;

  if IsProbablyPrimeInt(n) then return [[n],[]]; fi;
  k := 1;
  repeat k := 2*k; until f(k) >= n;
  step := k/4;
  while f(k) <> n and IsInt(step) do
    if f(k) > n then k := k - step; else k := k + step; fi;
    step := step/2;
  od;
  if f(k) <> n then return [[],[n]]; fi;
  val := List(DivisorsInt(k),f);
  fact := [n];
  for fk in val do
    for i in [1..Length(fact)] do
      gcd := Gcd(fact[i],fk);
      if not gcd in [1,fact[i]] then fact[i] := [gcd,fact[i]/gcd]; fi;
    od;
    fact := Flat(fact);
  od;
  return [Filtered(fact,IsProbablyPrimeInt),
          Filtered(fact,q->not IsProbablyPrimeInt(q))];
end;
MakeReadOnlyGlobal("FactorsMultFunc");


# Power Check

FactorsPowerCheck := function (n,SplittingFunction,SplittingFunctionName)

  local  m,k,FactorsOfm,factors;

  m := SmallestRootInt(n); 
  k := LogInt(n,m);
  if m < n 
  then Info(IntegerFactorizationInfo,1,n," = ",m,"^",k); fi;
  if IsProbablyPrimeInt(m) 
  then return [ListWithIdenticalEntries(k,m),[]];
  elif m = n then return [[],[n]]; 
  else 
    FactorsOfm := [[],[m]];
    ApplyFactoringMethod(SplittingFunction,[],FactorsOfm,infinity,
                         [SplittingFunctionName,
                          ", Number to be factored : ","n"]);
    factors := [Concatenation(ListWithIdenticalEntries(k,FactorsOfm[1])),
                Concatenation(ListWithIdenticalEntries(k,FactorsOfm[2]))];
    Sort(factors[1]); Sort(factors[2]);
    return factors;
  fi;
end;
MakeReadOnlyGlobal("FactorsPowerCheck");


# Check for factors VERY close to the square root of n

FactorsFermat := function ( n, steps )

  local  a, b, a2, b2, d,
         steps1, steps2, result;

  a := RootInt(n,2); a2 := a^2;
  d := 2*a + 1; b := 0;
  steps1 := 0; steps2 := 0;
  repeat
    if steps1 > steps then return [ [  ], [ n ] ]; fi;
    a      := a + 1;
    a2     := a2 + d;
    d      := d + 2;
    b2     := a2 - n;
    steps1 := steps1 + 1;
    if   not b2 mod 64 in [ 0, 1, 4, 9, 16, 17, 25, 33, 36, 41, 49, 57 ]
      or not b2 mod 45 in [ 0, 1, 4, 9, 10, 16, 19, 25, 31, 34, 36, 40 ]
      or not b2 mod  7 in [ 0, 1, 2, 4 ]
      or not b2 mod 11 in [ 0, 1, 3, 4, 5, 9 ]
    then continue; fi;
    b      := RootInt(b2,2);
    steps2 := steps2 + 1;
  until b^2 = b2;
  b := RootInt(b2,2);
  Info(InfoFactInt,2,"FactorsFermat: #Steps = ",steps1," / ", steps2);
  result := [];
  result[1] := Filtered([a-b,a+b],IsProbablyPrimeInt);
  result[2] := Difference([a-b,a+b],result[1]);
  return result;
end;
MakeReadOnlyGlobal("FactorsFermat");


# Check for n = b^k +/- 1

FactorsAurifeuillian := function ( n )

  local b, k, s, FactorsOfP, PolyFactors, factors, c, m, j, p, a;

  for c  in [ -1, 1 ]  do
    b := SmallestRootInt( n - c );
    if b < n - c  then
      k := LogInt( n - c, b );
      if c = -1  then
        s := " - 1";
      else
        s := " + 1";
      fi;
      Info( IntegerFactorizationInfo, 1, n, " = ", b, "^", k, s);
      if c = -1 then
        FactorsOfP := DivisorsInt(k);
      else
        FactorsOfP := Difference(DivisorsInt(2*k), DivisorsInt(k));
      fi;
      PolyFactors := List(FactorsOfP,
                          i -> ValuePol(CyclotomicPol(i), b));
      Info( IntegerFactorizationInfo, 1,
            "The factors corresponding to ", "polynomial factors are\n",
            PolyFactors );
      factors := [ [  ], [  ] ];
      for j in [1..Length(FactorsOfP)] do
        a := PolyFactors[j];
        if b <= Length(BRENTFACTORSAVAILABLE)
          and not IsBound(BRENTFACTORS[b]) and BRENTFACTORSAVAILABLE[b]
        then
          ReadPackage("factint",Concatenation("tables/brfac",String(b)));
        fi;
        if b <= Length(BRENTFACTORSAVAILABLE)
          and IsBound(BRENTFACTORS[b])
          and IsBound(BRENTFACTORS[b][FactorsOfP[j]])
        then
          for p in BRENTFACTORS[b][FactorsOfP[j]] do
            while a mod p = 0 do
              Add(factors[1],p);
              a := a/p;
            od;
          od;
        fi;
        if IsProbablyPrimeInt(a) then
          Add(factors[1], a);
        else
          Add(factors[2], a);
        fi;
      od;
      return factors;
    fi;
  od;
  return [ [  ], [ n ] ];
end;
MakeReadOnlyGlobal("FactorsAurifeuillian");


#############################################################################
##
#F  FactInt( <n> ) . . . . . . . . . . prime factorization of the integer <n>
#F                                                      (partial or complete)
##
##  Recognized options are:
##
##  <TDHints>          a list of additional trial divisors
##  <RhoSteps>         number of steps for Pollard's Rho
##  <RhoCluster>       interval for Gcd computation in Pollard's Rho
##  <Pminus1Limit1>    first stage limit for Pollard's $p-1$
##  <Pminus1Limit2>    second stage limit for Pollard's $p-1$ 
##  <Pplus1Residues>   number of residues to be tried in William's $p+1$
##  <Pplus1Limit1>     first stage limit for William's $p+1$
##  <Pplus1Limit2>     second stage limit for William's $p+1$
##  <ECMCurves>        number of elliptic curves to be tried by 
##                     the Elliptic Curves Method (ECM),
##                     also admissible: a function that takes the number to
##                     be factored and returns the desired number of curves 
##  <ECMLimit1>        initial first stage limit for ECM
##  <ECMLimit2>        initial second stage limit for ECM
##  <ECMDelta>         increment for first stage limit in ECM
##                     (the second stage limit is also incremented 
##                     appropriately)
##  <ECMDeterministic> if true, the choice of curves in ECM is deterministic,
##                     i.e. repeatable 
##  <FactIntPartial>   if true, the partial factorization obtained by
##                     applying the factoring methods whose time complexity 
##                     depends mainly on the size of the factors to be found
##                     and less on the size of <n> (see manual) is returned
##                     and the factor base methods (MPQS and CFRAC) are not
##                     used to complete the factorization for numbers that
##                     exceed the bound given by <CFRACLimit> resp.
##                     <MPQSLimit>; default: false
##  <FBMethod>         specifies which of the factor base methods should be
##                     used to do the ``hard work''; currently implemented:
##                     `"CFRAC"' and `"MPQS"'
##  <CFRACLimit>       specifies the maximal number of decimal digits of an
##                     integer to which the Continued Fraction Algorithm
##                     (CFRAC) should be applied (only used when 
##                     <FactIntPartial> is true)
##  <MPQSLimit>        as above, for the Multiple Polynomial Quadratic
##                     Sieve (MPQS)
##
##  `FactInt' returns a list of two lists, where the first list contains the
##  prime factors of <n> which have been found, and the second one contains
##  the remaining unfactored part(s), if there are any.
##
InstallGlobalFunction(FactInt,

function ( n )

  local  TDHints,RhoSteps,RhoCluster,
         Pminus1Limit1,Pminus1Limit2,
         Pplus1Residues,Pplus1Limit1,Pplus1Limit2,
         ECMCurves,ECMLimit1,ECMLimit2,ECMDelta,
         FactIntPartial,FBMethod,CFRACLimit,MPQSLimit,
         IsNonnegInt,StateInfo,LastMentioned,
         FactorizationObtainedSoFar,Result,sign,
         CFRACBound,MPQSBound,StartingTime,UsedTime,
         fib_res;

  IsNonnegInt := n->(IsInt(n) and n >= 0);

  StateInfo := function ()
    if FactorizationObtainedSoFar[2] <> [] 
      and not (    IsBound(LastMentioned) 
               and FactorizationObtainedSoFar[1] = LastMentioned) 
    then
      Info(IntegerFactorizationInfo,2,"");
      Info(IntegerFactorizationInfo,1,
           "Factors already found : ",FactorizationObtainedSoFar[1]);
      Info(IntegerFactorizationInfo,1,"");
      LastMentioned := ShallowCopy(FactorizationObtainedSoFar[1]);
    fi;
  end;

  if  not IsInt(n)
  then Error("Usage : FactInt( <n> ), for an integer <n>"); fi;

  if AbsInt(n) < 10^12 then
    Info(IntegerFactorizationInfo,3," | ",n," | ",
         "< 10^12, so use library function `FactorsInt'");
    return [FactorsInt(n),[]];
  fi;

  StartingTime := Runtime();

  # Get options / set default values

  TDHints := ValueOption("TDHints");
  if not IsList(TDHints) or not ForAll(TDHints,IsPosInt) 
  then TDHints := []; fi;
  RhoSteps := ValueOption("RhoSteps"); 
  if not IsNonnegInt(RhoSteps) then RhoSteps := 16384; fi;
  RhoCluster := ValueOption("RhoCluster"); 
  if not IsPosInt(RhoCluster) 
  then RhoCluster := Maximum(Minimum(Int(LogInt(AbsInt(n),2)^2/100),
                                     Int(RhoSteps/10)),16); fi;
  Pminus1Limit1 := ValueOption("Pminus1Limit1");
  if not IsNonnegInt(Pminus1Limit1) then Pminus1Limit1 := 10000; fi;
  Pminus1Limit2 := ValueOption("Pminus1Limit2");
  if not IsNonnegInt(Pminus1Limit2) 
  then Pminus1Limit2 := 40 * Pminus1Limit1; fi;
  Pplus1Residues := ValueOption("Pplus1Residues");
  if not IsNonnegInt(Pplus1Residues) then Pplus1Residues := 2; fi;
  Pplus1Limit1 := ValueOption("Pplus1Limit1");
  if not IsNonnegInt(Pplus1Limit1) then Pplus1Limit1 := 2000; fi;
  Pplus1Limit2 := ValueOption("Pplus1Limit2");
  if not IsNonnegInt(Pplus1Limit2) 
  then Pplus1Limit2 := 40 * Pplus1Limit1; fi;
  ECMCurves := ValueOption("ECMCurves");
  ECMLimit1 := ValueOption("ECMLimit1");
  ECMLimit2 := ValueOption("ECMLimit2");
  ECMDelta  := ValueOption("ECMDelta");
  FactIntPartial := ValueOption("FactIntPartial");
  if not IsBool(FactIntPartial) or FactIntPartial = fail 
  then FactIntPartial := false; fi;
  FBMethod := ValueOption("FBMethod");
  if not IsString(FBMethod) or not FBMethod in ["MPQS","CFRAC"]
  then FBMethod := "MPQS"; fi;
  CFRACLimit := ValueOption("CFRACLimit");
  if not IsPosInt(CFRACLimit) then CFRACLimit := 40; fi;
  MPQSLimit := ValueOption("MPQSLimit");
  if not IsPosInt(MPQSLimit) then MPQSLimit := 40; fi;

  if n < 0  then sign := -1; else sign := 1; fi;    

  FactorizationObtainedSoFar := [[],[AbsInt(n)]];

  # First of all, check whether n = b^k +- 1 for some b, k

  ApplyFactoringMethod(FactorsAurifeuillian,[],
                       FactorizationObtainedSoFar,infinity,
                       ["Check for n = b^k +/- 1"]);
  StateInfo();

  # Special case k! +/- 1

  if n mod 620448401733239439360000 in [1,620448401733239439359999] then
    if   IsEmpty(K_FACTORIAL_M1_FACTORS)
    then ReadPackage("factint","tables/factorial.g"); fi;
    if n mod 6 = 1 then
      ApplyFactoringMethod(FactorsTD,[K_FACTORIAL_P1_FACTORS],
                           FactorizationObtainedSoFar,infinity,
                           ["Trial division by factors of k!+1"]);
    else
      ApplyFactoringMethod(FactorsTD,[K_FACTORIAL_M1_FACTORS],
                           FactorizationObtainedSoFar,infinity,
                           ["Trial division by factors of k!-1"]);
    fi;
    StateInfo();
  fi;

  # Special case Primorial(k) +/- 1

  if n mod 32589158477190044730 in [1,32589158477190044729] then
    if   IsEmpty(K_PRIMORIAL_M1_FACTORS)
    then ReadPackage("factint","tables/primorial.g"); fi;
    if n mod 6 = 1 then
      ApplyFactoringMethod(FactorsTD,[K_PRIMORIAL_P1_FACTORS],
                           FactorizationObtainedSoFar,infinity,
                           ["Trial division by factors of Primorial(k)+1"]);
    else
      ApplyFactoringMethod(FactorsTD,[K_PRIMORIAL_M1_FACTORS],
                           FactorizationObtainedSoFar,infinity,
                           ["Trial division by factors of Primorial(k)-1"]);
    fi;
    StateInfo();
  fi;

  # Special case Fibonacci numbers

  fib_res := List([13,21,34,55,89,144], m -> n mod m);
  if ForAll([1..6],i->fib_res[i] in FIB_RES[i]) then
    ApplyFactoringMethod(FactorsMultFunc,[Fibonacci],
                         FactorizationObtainedSoFar,infinity,
                         ["Factors of Fibonacci(k) by divisors of k"]);
    if   IsEmpty(FACTORS_FIB)
    then ReadPackage("factint","tables/fibo.g"); fi;
    ApplyFactoringMethod(FactorsTD,[FACTORS_FIB],
                         FactorizationObtainedSoFar,infinity,
                         ["Trial division by factors of Fibonacci(k)"]);
    StateInfo();
  fi;

  # Special case 3^k - 2^k

  if n mod 2520 in [1,5,19,65,211,665,1051,1219,1265,1531,2059] then
    ApplyFactoringMethod(FactorsMultFunc,[k->3^k-2^k],
                         FactorizationObtainedSoFar,infinity,
                         ["Factors of 3^k-2^k by divisors of k"]);
    if   IsEmpty(POW3_M_POW2_FACTORS)
    then ReadPackage("factint","tables/3k2k.g"); fi;
    ApplyFactoringMethod(FactorsTD,[POW3_M_POW2_FACTORS],
                         FactorizationObtainedSoFar,infinity,
                         ["Trial division by factors of 3^k-2^k"]);
    StateInfo();
  fi;

  # The 'naive' methods

  ApplyFactoringMethod(FactorsTD,[],
                       FactorizationObtainedSoFar,infinity,
                       ["Trial division by all primes p < 1000"]);
  StateInfo();
  ApplyFactoringMethod(FactorsTD,[Primes2],
                       FactorizationObtainedSoFar,infinity,
                       ["Trial division by some already known primes"]);
  StateInfo();
  if TDHints <> [] then
  ApplyFactoringMethod(FactorsTD,[TDHints],
                       FactorizationObtainedSoFar,infinity,
                       ["Trial division by factors given as <TDHints>"]); fi;
  StateInfo();
  ApplyFactoringMethod(FactorsTD,
                       [Filtered(List(Filtered(Difference(Set(IDENTS_GVAR()),
                                               NAMES_SYSTEM_GVARS),ISB_GVAR),
                                      ValueGlobal),
                                 obj -> TNUM_OBJ_INT(obj)=1)],
                       FactorizationObtainedSoFar,infinity,
                       ["Trial division by user GVar's in workspace"]);
  StateInfo();
  ApplyFactoringMethod(FactorsPowerCheck,[FactInt,"FactInt"],
                       FactorizationObtainedSoFar,infinity,
                       ["Check for perfect powers"]);
  StateInfo();

  # Special case of two factors VERY close to the square root

  ApplyFactoringMethod(FactorsFermat,[100],
                       FactorizationObtainedSoFar,infinity,
                       ["Fermat's method"]);
  StateInfo();

  # Let 'FactorsRho', 'FactorsPminus1', 'FactorsPplus1' and 'FactorsECM' 
  # cast out the medium-sized factors

  if RhoSteps > 0 then
  ApplyFactoringMethod(FactorsRho,[1,RhoCluster,RhoSteps],
                       FactorizationObtainedSoFar,infinity,
                       ["Pollard's Rho\nSteps = ",RhoSteps,
                        ", Cluster = ",RhoCluster,
                        "\nNumber to be factored : ","n"]); fi;
  StateInfo();

  if Pminus1Limit1 > 0 then
  ApplyFactoringMethod(FactorsPminus1,[2,Pminus1Limit1,Pminus1Limit2],
                       FactorizationObtainedSoFar,infinity,
                       ["Pollard's p - 1\nLimit1 = ",
                        Pminus1Limit1,", Limit2 = ",Pminus1Limit2,
                        "\nNumber to be factored : ","n"]); fi;
  StateInfo();
  
  if Pplus1Residues > 0 and Pplus1Limit1 > 0 then
  ApplyFactoringMethod(FactorsPplus1,
                       [Pplus1Residues,Pplus1Limit1,Pplus1Limit2],
                       FactorizationObtainedSoFar,infinity,
                       ["Williams' p + 1\nResidues = ",Pplus1Residues,
                        ", Limit1 = ",Pplus1Limit1,", Limit2 = ",
                        Pplus1Limit2,
                        "\nNumber to be factored : ","n"]); fi;
  StateInfo();

  if ECMLimit1 > 0 and ECMCurves <> 0 then
  ApplyFactoringMethod(FactorsECM,[ECMCurves,ECMLimit1,ECMLimit2,ECMDelta],
                       FactorizationObtainedSoFar,infinity,
                       ["Elliptic Curves Method (ECM)\n",
                        "Curves = ",ECMCurves,"\nInit. Limit1 = ",
                        ECMLimit1,", Init. Limit2 = ",ECMLimit2,
                        ", Delta = ",ECMDelta,
                        "\nNumber to be factored : ","n"]); fi;
  StateInfo();

  # Let FactorsMPQS or FactorsCFRAC
  # do the really hard work, if <FactIntPartial> is false
  # or the remaining composite factors are smaller than
  # the upper bounds given by CFRACLimit and MPQSLimit 

  if not FactIntPartial 
  then CFRACBound := infinity;      MPQSBound := infinity;
  else CFRACBound := 10^CFRACLimit; MPQSBound := 10^MPQSLimit; fi;

  if FBMethod = "MPQS" then
  ApplyFactoringMethod(FactorsMPQS,[],
                       FactorizationObtainedSoFar,MPQSBound,
                       ["Multiple Polynomial Quadratic Sieve (MPQS)\n",
                       "Number to be factored : ","n"]:NoPreprocessing);
  elif FBMethod = "CFRAC" then
  ApplyFactoringMethod(FactorsCFRAC,[],
                       FactorizationObtainedSoFar,CFRACBound,
                       ["Continued Fraction Algorithm (CFRAC)\n",
                       "Number to be factored : ","n"]:NoPreprocessing);
  fi;

  Result := FactorizationObtainedSoFar;
  if Result[1] <> [] then Result[1][1] := Result[1][1] * sign;
                     else Result[2][1] := Result[2][1] * sign; fi;

  Info(IntegerFactorizationInfo,1,"");
  Info(IntegerFactorizationInfo,1,"The result is\n",Result,"\n");
  UsedTime := Runtime() - StartingTime;
  Info(IntegerFactorizationInfo,2,"The total runtime was ",
                                   TimeToString(UsedTime),"\n");
  FactorizationCheck(n,Result);
  return Result;
end);

#############################################################################
##
#F  IntegerFactorization( <n> ) . . . . . .  prime factors of the integer <n>
## 
##  Returns the list of prime factors of the integer <n>.
##
InstallGlobalFunction(IntegerFactorization,
function (n)
  if   not IsInt(n) 
  then Error("Usage : IntegerFactorization( <n> ), ",
             "where n has to be an integer"); fi;

  return FactInt(n:FactIntPartial:=false)[1];
end);

#############################################################################
##
#M  Factors( Integers, <n> )  . . . . . . . . . . factorization of an integer
##
InstallMethod( Factors,
               "FactInt: for integers", true, [ IsIntegers, IsInt ], 1,

  function ( Integers, n )
    if   IsSmallIntRep( n )
    then return FactorsInt( n );
    else return IntegerFactorization(n); fi;
  end );

#############################################################################
##
#E  general.gi . . . . . . . . . . . . . . . . . . . . . . . . . .  ends here