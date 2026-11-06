SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(4*10^9);tt:=Cputime();
load "work/deg10/continuation_20261029b/deps/pair_af_place_linear.m";
D:=eval Read("work/deg10/continuation_20261029b/C226/proven_picard.m");k:=GF(3);PP:=PolynomialRing(k);B<t>:=FunctionField(k);RF<Y>:=PolynomialRing(B);
function Rat(v) return B!(PP!v[1])/B!(PP!v[2]);end function;
f:=RF![Rat(v):v in D[1]];f/:=LeadingCoefficient(f);F<y>:=FunctionField(f);assert Genus(F) eq 25;
function GetFunction(data) return F![Rat(v):v in data];end function;
OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);fb:=[];
for j in [1..#D[4]] do OO:=D[2][j] eq 1 select OF else OI;I:=ideal<OO|[GetFunction(v):v in D[4][j]]>;pl:=Place(I);assert Degree(pl) eq D[3][j];Append(~fb,pl);end for;
Z:=FreeAbelianGroup(#D[5][1]);G,q:=quo<Z|[Z!v:v in D[6]]>;fc:=[q(Z!v):v in D[5]];A:=Divisor(fb[D[8]]);coords:=[GetFunction(v):v in D[7]];
function SparseIdeals(sp)
 IF:=ideal<OF|1>;II:=ideal<OI|1>;d:=0;
 for r in sp do j,e:=Explode(r);d+:=e*D[3][j];if D[2][j] eq 1 then IF*:=Ideal(fb[j])^e;else II*:=Ideal(fb[j])^e;end if;end for;return IF,II,d;
end function;
source:=eval Read("work/deg10/continuation_20261029b/C226/filter26.m");fields:=AssociativeArray();for r in source[4] do fields[r[1]]:=ext<k|PP!r[2]>;end for;
ww:=eval Read("work/deg10/continuation_20261029b/C226/all_picard_witnesses.m");out:=[];
for shard in [0..7] do
 data:=eval Read(Sprintf("work/deg10/continuation_20261029b/C226/picard%o.m",shard));wi:=0;smooth:=0;
 for im in data[2] do
  kk,cl:=Explode(im);d:=kk[1];ky:=kk[2..#kk];K:=fields[d];z:=K.1;bc:=[Evaluate(PP!ky[1..d],z),Evaluate(PP!ky[d+1..2*d],z)];pl:=C48PairAFPlaceLinear(coords,bc);assert Degree(pl) eq d;D0:=Divisor(pl)-d*A;j:=Index(fb,pl);
  if j gt 0 then cla:=fc[j]-d*fc[D[8]];smooth+:=1;
  else
   wi+:=1;w:=ww[shard+1][wi];assert w[2] eq 0;IE,JE,de:=SparseIdeals(w[4]);IR,JR,dr:=SparseIdeals(w[5]);IA,JA,da:=SparseIdeals(w[6]);assert de eq 25 and dr eq 25 and da eq 1;
   I0,J0:=Ideals(D0);AA,AB:=Ideals(A);assert IA eq AA and JA eq AB;h:=GetFunction(w[7]);assert h ne 0;assert IE/(IR*(h*OF)) eq I0 and JE/(JR*(h*OI)) eq J0;
   cla:=G!w[3];v:=G!0;for r in w[4] do v+:=r[2]*fc[r[1]];end for;for r in w[5] do v-:=r[2]*fc[r[1]];end for;assert cla eq v;
  end if;assert cla eq G!cl;
 end for;assert wi eq #ww[shard+1];Append(~out,[shard,#data[2],wi,smooth]);print "C226_PICARD_IDEAL_REPLAY",out[#out],Cputime(tt);
end for;
Write("work/deg10/continuation_20261029b/C226/native_ideal_replay.m",Sprint(out));print "C226_ALL_PICARD_IDEAL_REPLAYS_COMPLETE",Cputime(tt);quit;
