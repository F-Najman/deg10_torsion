SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(4*10^9);tt:=Cputime();
root:="work/deg10/C236_resolution_20261030/enumeration/";D:=eval Read(root cat "target18_proven_picard.m");k:=GF(5);PP:=PolynomialRing(k);B<t>:=FunctionField(k);RF<Y>:=PolynomialRing(B);
function Rat(v) return B!(PP!v[1])/B!(PP!v[2]);end function;
f:=RF![Rat(v):v in D[1]];f/:=LeadingCoefficient(f);F<y>:=FunctionField(f);assert Genus(F) eq 7;
function GetFunction(data) return F![Rat(v):v in data];end function;
OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);fb:=[];
for j in [1..#D[4]] do OO:=D[2][j] eq 1 select OF else OI;I:=ideal<OO|[GetFunction(v):v in D[4][j]]>;pl:=Place(I);assert Degree(pl) eq D[3][j];Append(~fb,pl);end for;
Z:=FreeAbelianGroup(#D[5][1]);G,q:=quo<Z|[Z!v:v in D[6]]>;fc:=[q(Z!v):v in D[5]];A:=Divisor(fb[D[8]]);AA,AB:=Ideals(A);
function SparseIdeals(sp)
 IF:=ideal<OF|1>;II:=ideal<OI|1>;d:=0;
 for r in sp do j,e:=Explode(r);d+:=e*D[3][j];if D[2][j] eq 1 then IF*:=Ideal(fb[j])^e;else II*:=Ideal(fb[j])^e;end if;end for;return IF,II,d;
end function;
data:=eval Read(root cat Sprintf("filter18_images_%o.m",TAG));ww:=eval Read(root cat Sprintf("filter18_%o_RR_witnesses.m",TAG));wi:=0;smooth:=0;
for im in data do
 kk,e,mult,cl,flag,ibasis,unweighted:=Explode(im);d:=kk[1];assert d eq e*mult;OO:=flag eq 1 select OF else OI;I:=ideal<OO|[GetFunction(v):v in ibasis]>;pl:=Place(I);assert Degree(pl) eq e;D0:=Divisor(pl)-e*A;j:=Index(fb,pl);
 if j gt 0 then cla:=fc[j]-e*fc[D[8]];smooth+:=1;
 else
  wi+:=1;w:=ww[wi];assert w[2] eq 0;IE,JE,de:=SparseIdeals(w[4]);IR,JR,dr:=SparseIdeals(w[5]);IA,JA,da:=SparseIdeals(w[6]);assert de eq 7 and dr eq 7 and da eq 1;
  I0,J0:=Ideals(D0);assert IA eq AA and JA eq AB;h:=GetFunction(w[7]);assert h ne 0;assert IE/(IR*(h*OF)) eq I0 and JE/(JR*(h*OI)) eq J0;
  cla:=q(Z!w[3]);v:=G!0;for r in w[4] do v+:=r[2]*fc[r[1]];end for;for r in w[5] do v-:=r[2]*fc[r[1]];end for;assert cla eq v;
 end if;assert cla eq q(Z!unweighted) and mult*cla eq q(Z!cl);
 if wi mod 1000 eq 0 and wi gt 0 then print "C236_IDEAL_REPLAY_PROGRESS",TAG,wi,Cputime(tt);end if;
end for;assert wi eq #ww;out:=[*TAG,#data,wi,smooth,Cputime(tt)*];Write(root cat Sprintf("native_ideal_replay_%o.m",TAG),Sprint(out):Overwrite:=true);print "C236_IDEAL_REPLAY_COMPLETE",out;quit;
