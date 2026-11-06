SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(4*10^9);tt:=Cputime();
root:="work/deg10/continuation_20261030b/C81/";D:=eval Read(root cat "x0_81_proven_picard.m");k:=GF(2);PP:=PolynomialRing(k);B<t>:=FunctionField(k);RF<Y>:=PolynomialRing(B);
function Rat(v) return B!(PP!v[1])/B!(PP!v[2]);end function;
f:=RF![Rat(v):v in D[1]];f/:=LeadingCoefficient(f);F<y>:=FunctionField(f);assert Genus(F) eq 4;
function GetFunction(data) return F![Rat(v):v in data];end function;
OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);fb:=[];
for j in [1..#D[4]] do OO:=D[2][j] eq 1 select OF else OI;I:=ideal<OO|[GetFunction(v):v in D[4][j]]>;pl:=Place(I);assert Degree(pl) eq D[3][j];Append(~fb,pl);end for;
Z:=FreeAbelianGroup(#D[5][1]);G,q:=quo<Z|[Z!v:v in D[6]]>;fc:=[q(Z!v):v in D[5]];A:=Divisor(fb[D[8]]);AA,AB:=Ideals(A);
function SparseIdeals(sp)
 IF:=ideal<OF|1>;II:=ideal<OI|1>;d:=0;
 for r in sp do j,e:=Explode(r);d+:=e*D[3][j];if D[2][j] eq 1 then IF*:=Ideal(fb[j])^e;else II*:=Ideal(fb[j])^e;end if;end for;return IF,II,d;
end function;

data:=eval Read(root cat "x0_filter.m");ww:=eval Read(root cat "x0_parent_witnesses.m");wi:=0;nr:=0;remain:=[];H:=sub<G|[q(Z!v):v in D[9]]>;assert Invariants(G) eq[3,18,0] and Invariants(H) eq[3,9];
for rec in data[2] do ci,d,ky,vals,rows,passed:=Explode(rec);assert #rows eq 1;anypass:=false;
 for im in rows do e,mult,cl,yes,flag,ibasis:=Explode(im);assert d eq e*mult and mult eq 1;OO:=flag eq 1 select OF else OI;I:=ideal<OO|[GetFunction(v):v in ibasis]>;assert IsPrime(I);pl:=Place(I);assert Degree(pl) eq e;D0:=Divisor(pl)-e*A;j:=Index(fb,pl);
 if j gt 0 then cla:=fc[j]-e*fc[D[8]];else
  wi+:=1;w:=ww[wi];assert w[2] eq 0;IE,JE,de:=SparseIdeals(w[4]);IR,JR,dr:=SparseIdeals(w[5]);IA,JA,da:=SparseIdeals(w[6]);assert de eq 4 and dr eq 4 and da eq 1;I0,J0:=Ideals(D0);assert IA eq AA and JA eq AB;h:=GetFunction(w[7]);assert h ne 0;assert IE/(IR*(h*OF)) eq I0 and JE/(JR*(h*OI)) eq J0;cla:=q(Z!w[3]);v:=G!0;for r in w[4] do v+:=r[2]*fc[r[1]];end for;for r in w[5] do v-:=r[2]*fc[r[1]];end for;assert cla eq v;
 end if;assert mult*cla eq q(Z!cl);assert yes eq(cla in H);anypass or:=yes;nr+:=1;
 end for;assert passed eq anypass;if passed then Append(~remain,ci);end if;
end for;assert wi eq #ww and nr eq 10 and remain eq[1,2,3,7,10];Write(root cat "parent_native_x0_filter.m",Sprint([*nr,wi,remain*]):Overwrite:=true);print "C81_PARENT_NATIVE_X0_FILTER_COMPLETE",nr,wi,remain,Cputime(tt);quit;
