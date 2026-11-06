SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(3*10^9);tt:=Cputime();
D:=eval Read("work/deg10/continuation_20261029b/C226/proven_picard.m");k:=GF(3);PP:=PolynomialRing(k);B<t>:=FunctionField(k);RF<Y>:=PolynomialRing(B);
function Rat(v) return B!(PP!v[1])/B!(PP!v[2]);end function;
f:=RF![Rat(v):v in D[1]];f/:=LeadingCoefficient(f);F<y>:=FunctionField(f);
function GetFunction(data) return F![Rat(v):v in data];end function;
OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);fb:=[];
for j in [1..#D[4]] do OO:=D[2][j] eq 1 select OF else OI;J:=ideal<OO|[GetFunction(v):v in D[4][j]]>;Append(~fb,Place(J));end for;
Z:=FreeAbelianGroup(#D[5][1]);G,q:=quo<Z|[Z!v:v in D[6]]>;fc:=[q(Z!v):v in D[5]];A:=D[8];ri:=D[9];b:=GetFunction(D[7][1]);labels:=[Valuation(b,fb[j]):j in ri];assert #Seqset(labels) eq 18;hgens:=[fc[j]-fc[A]:j in ri];Fr:=FreeAbelianGroup(18);h:=hom<Fr->G|hgens>;ker:=[Eltseq(Fr!v):v in Generators(Kernel(h))];
rows:=[* *];
for shard in [0..7] do
 ar:=eval Read(Sprintf("work/deg10/continuation_20261029b/C226/picard%o.m",shard));
 for r in ar[3] do if not r[5] or r[3] ne 0 then continue;end if;cl:=G!r[4];v:=Eltseq(cl@@h);assert h(Fr!v) eq cl;Append(~rows,[*r[1],r[2],r[4],v*]);end for;
end for;
Write("work/deg10/continuation_20261029b/C226/cusp_preimages.m",Sprint([*labels,Valuation(b,fb[A]),[Eltseq(v):v in hgens],D[6],ker,rows*]));
print "C226_GLOBAL_CUSP_PREIMAGES_COMPLETE",labels,Valuation(b,fb[A]),#rows,Cputime(tt);quit;
