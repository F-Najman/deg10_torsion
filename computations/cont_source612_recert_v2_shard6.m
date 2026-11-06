SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(8*10^9);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
Attach("work/deg10/continuation_20261029/deps/fast_hecke.m");
Attach("work/deg10/continuation_20261029/deps/fast_hecke_11.m");
Attach("work/deg10/continuation_20261029/deps/model_6_12.m");
function CuspBlocks(X)
C:=Curve(X);FF:=FunctionField(C);E:=EllipticCurve(X`_E);Q:=E![0,0];px:=FF!X`_P[1];ai:=aInvariants(E);R<Y>:=PolynomialRing(FF);
ys:=Roots(Y^2+(ai[1]*px+ai[3])*Y-px^3-ai[2]*px^2-ai[4]*px-ai[5]);assert #ys eq 2;P:=E![px,ys[1][1]];
z:=PrimitiveElement(BaseRing(X))^((#BaseRing(X)-1) div 6);if WeilPairing(P,2*Q,6) ne z then P:=-P;end if;assert WeilPairing(P,2*Q,6) eq z;
fs:=[-ai[2],1-ai[1],Discriminant(E)];anchor:=P[1];
for i in [0..5] do for j in [0..11] do
 A:=i*P+j*Q;if A eq E!0 then continue;end if;
 if A ne Q and A ne -Q then assert A[1] ne 0;Append(~fs,A[1]);end if;
 if A ne P and A ne -P then assert A[1] ne anchor;Append(~fs,A[1]-anchor);end if;
end for;end for;
cs:=Cusps(X);groups:=AssociativeArray();
for c in cs do
 sig:=[Valuation(f,c):f in fs];if IsDefined(groups,sig) then Append(~groups[sig],c);else groups[sig]:=[c];end if;
end for;
sigs:=Sort(Setseq(Keys(groups)));blocks:=[&+[Divisor(c):c in groups[sig]]:sig in sigs];ds:=[Degree(b):b in blocks];assert #blocks eq 24;
return blocks,sigs,ds;
end function;
load "work/deg10/continuation_20261029/deps/class_group_direct_RR_m3_18_root.m";
a:=eval Read("work/deg10/continuation_20261029/diagonal612_D5_source_checkpoint_v5.m");
X:=MDX612(GF(7));C:=Curve(X);blocks,sigs,ds:=CuspBlocks(X);assert sigs eq a[1] and ds eq a[2];c0:=Support(blocks[Index(ds,1)])[1];AF:=AlgorithmicFunctionField(FunctionField(C));
G,_,native:=ClassGroup(AF:ReductionDivisor:=FunctionFieldDivisor(Divisor(c0)),Proof:=true);
GG,_,ph:=MDClassGroupDirectRRCheckedV3(AF,"work/deg10/continuation_20261029/cont_source612_recert_v2_shard6_identities.log");assert GG eq G;
m2:=map<DivisorGroup(C)->G|D:->ph(FunctionFieldDivisor(D))>;cg:=[m2(blocks[i]-ds[i]*Divisor(c0)):i in [1..24]];
OB:=FreeAbelianGroup(#a[5][1][1]);OG,op:=quo<OB|[OB!r:r in a[4][1]]>;oc:=[op(OB!v):v in a[5][1]];point:=AssociativeArray();
for row in a[8] do
 key:=[row[1],row[2]];cl:=op(OB!row[4])-&+[row[3][j]*oc[j]:j in [1..24]];
 if IsDefined(point,key) then assert point[key] eq cl;else point[key]:=cl;end if;
end for;assert #Keys(point) eq 3828;
xs:=NoncuspidalPlaces(X,3);oe:=[point[[3,i]]:i in [1,2]];nextra:=[m2(Divisor(xs[i])-3*Divisor(c0)):i in [1,2]];
F:=FreeAbelianGroup(26);oh:=hom<F->OG|oc cat oe>;nh:=hom<F->G|cg cat nextra>;assert #Image(oh) eq 95551488 and #Image(nh) eq 95551488;
for v in Generators(Kernel(oh)) do assert nh(F!v) eq G!0;end for;for v in Generators(Kernel(nh)) do assert oh(F!v) eq OG!0;end for;
AB:=FreeAbelianGroup(#a[6]);ha:=hom<AB->OG|[op(OB!v):v in a[6]]>;ia:=hom<AB->OG|[op(OB!v):v in a[7]]>;
assert #Image(ha) eq 95551488;for v in Generators(Kernel(ha)) do assert ia(AB!v) eq OG!0;end for;
for c in oc do assert ia(c@@ha) eq OG!0;end for;
for i in [1,2] do
 x:=xs[i];A:=HeckeOperator612(X,13,x)-Divisor(x)-13*Divisor(Diamond612(X,13,x));assert Degree(A) eq 0;assert m2(A) eq nh(F!(ia(oe[i]@@ha)@@oh));
end for;
records:=[* *];
for d in [3..5] do
 xs:=NoncuspidalPlaces(X,d);assert #xs eq [96,468,3264][d-2];
 for i in [1..#xs] do
  if (i-1) mod 8 ne 6 then continue;end if;
  actual:=m2(Divisor(xs[i])-d*Divisor(c0));expected:=nh(F!(point[[d,i]]@@oh));assert actual eq expected;Append(~records,[*d,i,Eltseq(OB!(point[[d,i]]@@op))*]);
 end for;
 print "DECIC612_SOURCE_RECERT_PROGRESS",d,#records,Cputime(tt);
end for;
Write("work/deg10/continuation_20261029/cont_source612_recert_v2_shard6_arithmetic.txt",Sprint(records));print "DECIC612_SOURCE_RECERT_COMPLETE",6,#records,Cputime(tt);quit;
