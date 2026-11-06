SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(6*10^9);tt:=Cputime();AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");load "work/deg10/continuation_20261031a/C70/inputs/pair_af_place_linear.m";
root:="work/deg10/continuation_20261031a/C70/";k:=GF(3);PK:=PolynomialRing(k);X:=MDX1(35,k);F0,mp:=AlgorithmicFunctionField(FunctionField(Curve(X)));F:=RationalExtensionRepresentation(F0);assert Genus(F)eq 25;coords:=[F!mp(v):v in X`_coordinates_xy];OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);print "D31A_C70_TARGET35_READY",Degree(F),Genus(F),Cputime(tt);
function Get(E,Q)return C48PairAFPlaceLinear(coords,TateNormalForm_xy(E,Q));end function;
function ER(v)return [[Integers()!a:a in Coefficients(Numerator(v))],[Integers()!a:a in Coefficients(Denominator(v))]];end function;
function EF(v)return [ER(a):a in Eltseq(F!v)];end function;
function PD(p)return [*Degree(p),Order(Ideal(p))eq OF select 1 else 0,[EF(v):v in Basis(Ideal(p))]*];end function;
function Test(E,Q,d)
 assert Order(Q)eq 35;ps:=Get(E,Q);pd:=Get(E,2*Q);assert d mod Degree(ps)eq 0 and 2*d mod Degree(pd)eq 0;ns:=d div Degree(ps);nd:=2*d div Degree(pd);D:=-ns*Divisor(ps)-nd*Divisor(pd);K:=BaseRing(E);RX<xx>:=PolynomialRing(K);psi:=DivisionPolynomial(E,2);assert Degree(psi)eq 3 and IsSquarefree(psi);fac:=Factorization(psi);rows:=[* *];
 for ff in fac do po:=ff[1];e:=Degree(po);if e eq 1 then L:=K;s:=-Coefficient(po,0)/Coefficient(po,1);else L:=ext<K|po>;s:=L.1;end if;EL:=ChangeRing(E,L);QL:=EL![L!Q[1],L!Q[2]];RL<zz>:=PolynomialRing(L);EG,phi:=IsogenyFromKernel(EL,zz-s);QT:=phi(QL);assert Order(QT)eq 35;pt:=Get(EG,QT);assert d*e mod Degree(pt)eq 0;n:=d*e div Degree(pt);D+:=n*Divisor(pt);Append(~rows,[*e,n,PD(pt)*]);end for;assert &+[r[1]:r in rows]eq 3 and Degree(D)eq 0;h0:=Dimension(D);return [*PD(ps),ns,PD(pd),nd,rows,h0*];
end function;
src:=eval Read(root cat "source70_p3.m");fields:=AssociativeArray();for row in src[1]do fields[row[1]]:=ext<k|PK!row[2]>;end for;tests:=[* *];done:=0;
for ci in CASES do d,ky:=Explode(src[3][ci]);K:=fields[d];a:=K.1;b:=Evaluate(PK!ky[1..d],a);c:=Evaluate(PK!ky[d+1..2*d],a);E:=EllipticCurve([K|1-c,-b,-b,0,0]);P:=E![0,0];assert Order(P)eq 70;Q:=2*P;v:=Test(E,Q,d);mps:=[*v*];
 if v[6]ne 0 then T:=35*P;RX<xx>:=PolynomialRing(K);EE,phi:=IsogenyFromKernel(E,xx-T[1]);QQ:=phi(Q);Append(~mps,Test(EE,QQ,d));end if;
 Append(~tests,[*ci,d,mps*]);done+:=1;if done mod 10 eq 0 then Write(root cat TAG cat "_progress.m",Sprint(tests):Overwrite:=true);print "D31A_C70_A2_PROGRESS",TAG,done,#[v:v in tests|&and[r[6]ne 0:r in v[3]]],Cputime(tt);end if;
end for;
out:=[*[ER(v):v in Coefficients(DefiningPolynomial(F))],[EF(v):v in coords],tests*];Write(root cat TAG cat ".m",Sprint(out):Overwrite:=true);print "D31A_C70_A2_COMPLETE",TAG,#tests,[[v[1],v[2]]:v in tests|&and[r[6]ne 0:r in v[3]]],Cputime(tt);quit;
