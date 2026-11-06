SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(4*10^9);tt:=Cputime();
load "work/deg10/C232_resolution_20261030/deps/pair_af_place_linear.m";
a:=eval Read("work/deg10/C232_resolution_20261030/target_proven_picard.m");k:=GF(3);PP:=PolynomialRing(k);B<t>:=FunctionField(k);RF<Y>:=PolynomialRing(B);
function Rat(v) return B!(PP!v[1])/B!(PP!v[2]);end function;
f:=RF![Rat(v):v in a[1]];f/:=LeadingCoefficient(f);F<y>:=FunctionField(f);assert Genus(F) eq 17;
function GetFunction(v) return F![Rat(w):w in v];end function;
coords:=[GetFunction(v):v in a[7]];OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);fb:=[];
for j in [1..#a[4]] do OO:=a[2][j] eq 1 select OF else OI;J:=ideal<OO|[GetFunction(v):v in a[4][j]]>;pl:=Place(J);assert Degree(pl) eq a[3][j];Append(~fb,pl);end for;
Z:=FreeAbelianGroup(#a[6][1]);G,q:=quo<Z|[Z!v:v in a[6]]>;fc:=[q(Z!v):v in a[5]];A:=Divisor(fb[a[8]]);
function SparseIdeals(sp)
 IF:=ideal<OF|1>;II:=ideal<OI|1>;d:=0;
 for r in sp do j,e:=Explode(r);d+:=e*a[3][j];if a[2][j] eq 1 then IF*:=Ideal(fb[j])^e;else II*:=Ideal(fb[j])^e;end if;end for;return IF,II,d;
end function;
images:=(eval Read(Sprintf("work/deg10/C232_resolution_20261030/filter%o.m",TAG)))[1];fields:=AssociativeArray();for d in [6..10] do src:=eval Read(Sprintf("work/deg10/C232_resolution_20261030/source%o.m",d));fields[d]:=ext<k|PP!src[2]>;end for;
wi:=0;smooth:=0;
for im in images do
 ky:=im[1];d:=ky[1];e:=im[2];norm:=im[3];K:=fields[d];z:=K.1;xy:=[Evaluate(PP!ky[2..d+1],z),Evaluate(PP!ky[d+2..2*d+1],z)];pl:=C48PairAFPlaceLinear(coords,xy);assert Degree(pl) eq e and e*norm eq d;OO:=im[5] eq 1 select OF else OI;J:=ideal<OO|[GetFunction(v):v in im[6]]>;assert Ideal(pl) eq J;D0:=Divisor(pl)-e*A;j:=Index(fb,pl);
 if j gt 0 then cla:=fc[j]-e*fc[a[8]];smooth+:=1;
 else
  wi+:=1;w:=ReplayWitnesses[wi];assert w[2] eq 0;IE,JE,de:=SparseIdeals(w[4]);IR,JR,dr:=SparseIdeals(w[5]);IA,JA,da:=SparseIdeals(w[6]);assert de eq 17 and dr eq 17 and da eq 1;I0,J0:=Ideals(D0);AA,AB:=Ideals(A);assert IA eq AA and JA eq AB;h:=GetFunction(w[7]);assert h ne 0;assert IE/(IR*(h*OF)) eq I0 and JE/(JR*(h*OI)) eq J0;
  cla:=q(Z!w[3]);v:=G!0;for r in w[4] do v+:=r[2]*fc[r[1]];end for;for r in w[5] do v-:=r[2]*fc[r[1]];end for;assert cla eq v;
 end if;
 assert cla eq q(Z!im[7]) and norm*cla eq q(Z!im[4]);
end for;
assert wi eq #ReplayWitnesses;Write(Sprintf("work/deg10/C232_resolution_20261030/ideal_replay%o.m",TAG),Sprint([*TAG,#images,wi,smooth*]):Overwrite:=true);print "C232_IDEAL_REPLAY_COMPLETE",TAG,#images,wi,smooth,Cputime(tt);quit;
