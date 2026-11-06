SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(4*10^9);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
load "work/deg10/continuation_20261029b/deps/pair_af_place_linear.m";
load "work/deg10/continuation_20261029b/C240/portable/portable_inputs.m";
k:=GF(3);PP:=PolynomialRing(k);B<t>:=FunctionField(k);RF<Y>:=PolynomialRing(B);
function Rat(v) return B!(PP!v[1])/B!(PP!v[2]);end function;
f:=RF![Rat(v):v in PortableFieldEquation];f/:=LeadingCoefficient(f);F<y>:=FunctionField(f);assert Genus(F) eq 25;
function GetFunction(data) return F![Rat(v):v in data];end function;
coords:=[GetFunction(v):v in PortableSutherland];OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);fb:=[];
for j in [1..#PortableIdealBases] do OO:=PortableFiniteFlags[j] eq 1 select OF else OI;J:=ideal<OO|[GetFunction(v):v in PortableIdealBases[j]]>;pl:=Place(J);assert Degree(pl) eq PortableDegrees[j];Append(~fb,pl);end for;
Z:=FreeAbelianGroup(9);G,q:=quo<Z|[Z!v:v in PortableRelations]>;fc:=[q(Z!v):v in PortableClasses];A:=Divisor(fb[609]);
function SparseIdeals(sp)
 IF:=ideal<OF|1>;II:=ideal<OI|1>;d:=0;
 for r in sp do j,e:=Explode(r);d+:=e*PortableDegrees[j];if PortableFiniteFlags[j] eq 1 then IF*:=Ideal(fb[j])^e;else II*:=Ideal(fb[j])^e;end if;end for;
 return IF,II,d;
end function;
data:=eval Read(Sprintf("work/deg10/continuation_20261029b/C240/portable/shard%o_arithmetic.txt",SHARD));images:=data[2];source:=eval Read("work/deg10/continuation_20261029b/C240/filter_quotient20_linear_arithmetic.txt");fields:=AssociativeArray();
for r in source[2] do fields[r[1]]:=ext<k|PP!r[2]>;end for;
wi:=0;smooth:=0;
for im in images do
 kk,e,norm,cl:=Explode(im);mapno:=kk[1];d:=kk[2];ky:=kk[3..#kk];assert mapno eq 1;K:=fields[d];z:=K.1;b:=Evaluate(PP!ky[1..d],z);c:=Evaluate(PP!ky[d+1..2*d],z);E:=EllipticCurve([K|0,c,0,b*(1-b-c),0]);Q:=E![b,b];assert 40*Q eq E!0 and 20*Q ne E!0 and 8*Q ne E!0;
 xy:=TateNormalForm_xy(E,Q);pl:=C48PairAFPlaceLinear(coords,xy);assert Degree(pl) eq e and norm*e eq d;D0:=Divisor(pl)-e*A;j:=Index(fb,pl);
 if j gt 0 then cla:=fc[j]-e*fc[609];smooth+:=1;
 else
  wi+:=1;w:=ReplayWitnesses[wi];assert w[2] eq Degree(D0) and w[2] eq 0;
  IE,JE,de:=SparseIdeals(w[4]);IR,JR,dr:=SparseIdeals(w[5]);IA,JA,da:=SparseIdeals(w[6]);assert de eq 25 and dr eq 25 and da eq 1;
  I0,J0:=Ideals(D0);AA,AB:=Ideals(A);assert IA eq AA and JA eq AB;h:=GetFunction(w[7]);assert h ne 0;
  assert IE/(IR*(h*OF)) eq I0 and JE/(JR*(h*OI)) eq J0;
  cla:=G!w[3];v:=G!0;for r in w[4] do v+:=r[2]*fc[r[1]];end for;for r in w[5] do v-:=r[2]*fc[r[1]];end for;assert cla eq v;
 end if;
 assert norm*cla eq G!cl;
end for;
assert wi eq #ReplayWitnesses;
Write(Sprintf("work/deg10/continuation_20261029b/C240/portable/replay%o_result.txt",SHARD),Sprint([SHARD,#images,wi,smooth]));
print "PICARD_IDEAL_REPLAY_COMPLETE",SHARD,#images,wi,smooth,Cputime(tt);quit;
