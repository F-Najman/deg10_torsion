SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(4*10^9);tt:=Cputime();
counts:=[];records:=[* *];fielddefs:=[* *];curveStats:=[* *];
function Key(j)
 return Minimum([Eltseq(j^(2^i),GF(2)):i in [0..Degree(Parent(j),GF(2))-1]]);
end function;
for d in [1..10] do
 F<a>:=GF(2^d);look:=AssociativeArray();qual:=0;groups:={@ @};
 for j in F do
  for E in Twists(EllipticCurveFromjInvariant(j)) do
   if #E mod 77 ne 0 then continue;end if;
   assert j ne 0;qual+:=1;A,mp:=AbelianGroup(E);iv:=Invariants(A);Include(~groups,iv);assert &*[GCD(n,77):n in iv] eq 77;
   idx:=[i:i in [1..#iv]|iv[i] mod 77 eq 0][1];P:= (iv[idx] div 77)*mp(A.idx);assert Order(P) eq 77;
   k:=Key(j);if IsDefined(look,k) then continue;end if;
   assert Degree(MinimalPolynomial(j,GF(2))) eq d;
   look[k]:=[*d,[Integers()!c:c in Coefficients(DefiningPolynomial(F,GF(2)))],[[Integers()!c:c in Eltseq(v,GF(2))]:v in aInvariants(E)],[[Integers()!c:c in Eltseq(P[i],GF(2))]:i in [1..2]],k,#E,iv*];
  end for;
 end for;
 rat:=30*qual;prior:=&+[Integers()|e*counts[e]:e in [1..d-1]|d mod e eq 0];assert (rat-prior) mod d eq 0;Append(~counts,(rat-prior) div d);
 for k in Sort(Setseq(Keys(look))) do Append(~records,look[k]);end for;
 Append(~curveStats,[*d,qual,#Keys(look),[g:g in groups]*]);
 printf "C77_CENSUS_DEGREE d=%o qualifying=%o jorbits=%o NCclosed=%o groups=%o cpu=%o\n",d,qual,#Keys(look),counts[d],groups,Cputime(tt);
end for;
assert counts[1..9] eq [0,0,0,0,0,0,0,0,0];
Write("work/deg10/continuation_20261030b/C77/census.m",Sprint([*counts,curveStats,records*]):Overwrite:=true);
printf "C77_CENSUS_COMPLETE counts=%o records=%o cpu=%o\n",counts,#records,Cputime(tt);quit;
