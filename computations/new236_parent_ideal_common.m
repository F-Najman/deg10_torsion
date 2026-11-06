SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(4*10^9);tt:=Cputime();
load "work/deg10/C236_resolution_20261030/geometry/deps/pair_af_place_linear.m";
load "work/deg10/C236_resolution_20261030/geometry/quotient9_projection_helpers.m";
a:=eval Read(SNAPSHOT);k:=GF(5);PP:=PolynomialRing(k);B<t>:=FunctionField(k);RF<Y>:=PolynomialRing(B);
function Rat(v) return B!(PP!v[1])/B!(PP!v[2]);end function;
f:=RF![Rat(v):v in a[1]];f/:=LeadingCoefficient(f);F<y>:=FunctionField(f);assert Genus(F) eq GENUS;
function GetFunction(v) return F![Rat(w):w in v];end function;
coords:=[GetFunction(v):v in a[7]];OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);fb:=[];
for j in [1..#a[4]] do OO:=a[2][j] eq 1 select OF else OI;J:=ideal<OO|[GetFunction(v):v in a[4][j]]>;pl:=Place(J);assert Degree(pl) eq a[3][j];Append(~fb,pl);end for;
Z:=FreeAbelianGroup(#a[6][1]);G,q:=quo<Z|[Z!v:v in a[6]]>;fc:=[q(Z!v):v in a[5]];A0:=Divisor(fb[a[8]]);
OO:=a[12][1] eq 1 select OF else OI;Ap:=Place(ideal<OO|[GetFunction(v):v in a[12][3]]>);assert Degree(Ap) eq 1;A:=Divisor(Ap);ja:=Index(fb,Ap);assert ja gt 0;
function SparseIdeals(sp)
 IF:=ideal<OF|1>;II:=ideal<OI|1>;d:=0;
 for r in sp do j,e:=Explode(r);d+:=e*a[3][j];if a[2][j] eq 1 then IF*:=Ideal(fb[j])^e;else II*:=Ideal(fb[j])^e;end if;end for;return IF,II,d;
end function;
data:=eval Read(INPUT);images:=data[2];fields:=AssociativeArray();for row in data[1] do fields[row[1]]:=ext<k|PP!row[2]>;end for;
wi:=0;smooth:=0;idx:=0;
for im in images do
 idx+:=1;
 if GENUS eq 9 then d:=im[1];e:=im[3];mult:=im[4];ans:=im[5];K:=fields[d];z:=K.1;vals:=[Evaluate(PP!im[6][i*d+1..(i+1)*d],z):i in [0..2]];pl:=InverseThree(coords,vals);
 else d:=im[2];e:=im[4];mult:=im[5];ans:=im[6];K:=fields[d];z:=K.1;vals:=[Evaluate(PP!v,z):v in im[7]];pl:=C48PairAFPlaceLinear(coords,vals);end if;
 assert Degree(pl) eq e and e*mult eq d;D0:=Divisor(pl)-e*A;j:=Index(fb,pl);
 if j gt 0 then cla:=fc[j]-e*fc[ja];smooth+:=1;
 else
  wi+:=1;w:=ReplayWitnesses[wi];assert w[2] eq 0;IE,JE,de:=SparseIdeals(w[4]);IR,JR,dr:=SparseIdeals(w[5]);IA,JA,da:=SparseIdeals(w[6]);assert de eq GENUS and dr eq GENUS and da eq 1;I0,J0:=Ideals(D0);AA,AB:=Ideals(A0);assert IA eq AA and JA eq AB;h:=GetFunction(w[7]);assert h ne 0;assert IE/(IR*(h*OF)) eq I0 and JE/(JR*(h*OI)) eq J0;
  cla:=q(Z!w[3]);v:=G!0;for r in w[4] do v+:=r[2]*fc[r[1]];end for;for r in w[5] do v-:=r[2]*fc[r[1]];end for;assert cla eq v;
 end if;
 assert mult*cla eq q(Z!ans);
 if idx mod 100 eq 0 then print "C236_NATIVE_IDEAL_REPLAY_PROGRESS",TAG,idx,wi,smooth,Cputime(tt);end if;
end for;
assert wi eq #ReplayWitnesses;Write(OUTPUT,Sprint([*TAG,#images,wi,smooth*]):Overwrite:=true);print "C236_NATIVE_IDEAL_REPLAY_COMPLETE",TAG,#images,wi,smooth,Cputime(tt);quit;
