SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(3*10^9);tt:=Cputime();
D:=eval Read("work/deg10/continuation_20261029b/C226/proven_picard.m");k:=GF(3);PP:=PolynomialRing(k);B<t>:=FunctionField(k);RF<Y>:=PolynomialRing(B);
function Rat(v) return B!(PP!v[1])/B!(PP!v[2]);end function;
f:=RF![Rat(v):v in D[1]];f/:=LeadingCoefficient(f);F<y>:=FunctionField(f);
function GetFunction(data) return F![Rat(v):v in data];end function;
OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);fb:=[];
for j in [1..#D[4]] do OO:=D[2][j] eq 1 select OF else OI;J:=ideal<OO|[GetFunction(v):v in D[4][j]]>;Append(~fb,Place(J));end for;
Z:=FreeAbelianGroup(#D[5][1]);G,q:=quo<Z|[Z!v:v in D[6]]>;fc:=[q(Z!v):v in D[5]];A:=D[8];ri:=D[9];b:=GetFunction(D[7][1]);labels:=[Valuation(b,fb[j]):j in ri];assert #Seqset(labels) eq 18;hgens:=[fc[j]-fc[A]:j in ri];Fr:=FreeAbelianGroup(18);h:=hom<Fr->G|hgens>;ker:=[Eltseq(Fr!v):v in Generators(Kernel(h))];

nonrat:=[* *];for j in D[10] do Append(~nonrat,[*Valuation(b,fb[j]),Eltseq(fc[j]-3*fc[A])*]);end for;
vals:=Setseq({r[1]:r in nonrat});assert #vals eq 3;glob:=[* *];
for v in vals do rows:=[r:r in nonrat|r[1] eq v];assert #rows eq 2;cl:=G!rows[1][2]+G!rows[2][2];Append(~glob,[*v,Eltseq(cl)*]);end for;
Write("work/deg10/C226_resolution_20261030/all_cusp_classes.m",Sprint([*labels,[Eltseq(v):v in hgens],glob,D[6]*]));print "Q226_GLOBAL_CUSP_ORBIT_CLASSES_COMPLETE",glob,Cputime(tt);quit;
