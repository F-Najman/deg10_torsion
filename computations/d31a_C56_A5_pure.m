SetColumns(0);SetNthreads(1);SetQuitOnError(true);SetMemoryLimit(9*10^9);Alarm(1000);tt:=Cputime();AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");Attach("work/deg10/continuation_20261031a/C56/inputs/fast_hecke.m");load "work/deg10/continuation_20261031a/C56/inputs/pair_af_place_linear.m";
root:="work/deg10/continuation_20261031a/C56/";k:=GF(3);PK:=PolynomialRing(k);X:=MDX1(56,k);F0,mp:=AlgorithmicFunctionField(FunctionField(Curve(X)));F:=RationalExtensionRepresentation(F0);assert Genus(F)eq 61;coords:=[F!mp(v):v in X`_coordinates_xy];OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);print "D31A_C56_FULL_SOURCE_READY",Degree(F),Genus(F),Cputime(tt);
function Get(E,Q)return C48PairAFPlaceLinear(coords,TateNormalForm_xy(E,Q));end function;
function ER(v)return [[Integers()!a:a in Coefficients(Numerator(v))],[Integers()!a:a in Coefficients(Denominator(v))]];end function;
function EF(v)return [ER(a):a in Eltseq(F!v)];end function;
function PD(p)return [*Degree(p),Order(Ideal(p))eq OF select 1 else 0,[EF(v):v in Basis(Ideal(p))]*];end function;
function Make(E,Q,d)
 assert 56*Q eq E!0 and 28*Q ne E!0 and 8*Q ne E!0;p:=Get(E,Q);pd:=Get(E,5*Q);assert d mod Degree(p)eq 0 and 5*d mod Degree(pd)eq 0;ns:=d div Degree(p);nd:=5*d div Degree(pd);D:=-ns*Divisor(p)-nd*Divisor(pd);rows:=[* *];
 for z in MDIsogeniesFast(E,Q,5)do EE,QQ,e:=Explode(z);pt:=Get(EE,QQ);assert d*e mod Degree(pt)eq 0;n:=d*e div Degree(pt);D+:=n*Divisor(pt);Append(~rows,[*e,n,PD(pt)*]);end for;assert Degree(D)eq 0;return D,[*PD(p),ns,PD(pd),nd,rows*];
end function;
a:=eval Read(root cat "fourmap_filter.m");fs:=AssociativeArray();for z in a[1]do fs[z[1]]:=ext<k|PK!z[2]>;end for;cache:=AssociativeArray();obs:=AssociativeArray();out:=[* *];left:=[* *];
Write(root cat "A5_full_source_model.m",Sprint([*[ER(v):v in Coefficients(DefiningPolynomial(F))],[EF(v):v in coords]*]):Overwrite:=true);
for z in a[4]do if z[3]ne 0 then continue;end if;D:=DivisorGroup(F)!0;rs:=[* *];for si in z[2]do if not IsDefined(cache,si)then r:=a[2][si];d:=r[2];ky:=r[3];K:=fs[d];b:=Evaluate(PK!ky[1..d],K.1);c:=Evaluate(PK!ky[d+1..2*d],K.1);E:=EllipticCurve([K|1-c,-b,-b,0,0]);Q:=E![0,0];DD,oo:=Make(E,Q,d);cache[si]:=DD;obs[si]:=oo;end if;D+:=cache[si];Append(~rs,[*si,obs[si]*]);end for;h:=Dimension(D);Append(~out,[*z[1],z[2],rs,h*]);if h ne 0 then Append(~left,z);end if;
 Write(root cat "A5_pure_progress.m",Sprint([*out,left*]):Overwrite:=true);print "D31A_C56_A5_PURE_CASE",z[1],h,#out,#left,Cputime(tt);
end for;
Write(root cat "A5_pure.m",Sprint([*out,left*]):Overwrite:=true);print "D31A_C56_A5_PURE_COMPLETE",#out,#left,Cputime(tt);quit;
