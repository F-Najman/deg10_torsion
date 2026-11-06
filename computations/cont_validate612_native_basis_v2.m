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
GG,_,ph:=MDClassGroupDirectRRCheckedV3(AF,"work/deg10/continuation_20261029/612_native_basis_validation_v2_initial_identities.log");assert GG eq G;
m2:=map<DivisorGroup(C)->G|D:->ph(FunctionFieldDivisor(D))>;cg:=[m2(blocks[i]-ds[i]*Divisor(c0)):i in [1..24]];
F0:=RationalExtensionRepresentation(AF);GC,_,fwdmap:=MDClassGroupDirectRRCheckedV3(F0,"work/deg10/continuation_20261029/612_native_basis_validation_v2_identities.log");assert GC eq G;
fb:=[v`place:v in F0`cld`fb];high:=F0`cld`ClToDivHigh;checks:=[* *];
for i in [1..#fb] do
 D:=Divisor(fb[i]);cl:=fwdmap(D);T:=high(cl);delta:=D-T;assert Degree(delta) eq 0;ok,f:=IsPrincipal(delta);assert ok and Divisor(f) eq delta;
 Append(~checks,[*i,Degree(fb[i]),Eltseq(cl),true*]);
end for;
Write("work/deg10/continuation_20261029/612_native_basis_validation_v2.m",Sprint(checks));print "DECIC612_NATIVE_BASIS_VALIDATION_COMPLETE",#checks,Cputime(tt);quit;
