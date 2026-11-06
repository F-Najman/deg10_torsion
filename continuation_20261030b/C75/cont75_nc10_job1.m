JobIndex:=1;CaseIndices:=[1, 2];
/* New canonical degree20 quotient model attempt; no exclusion claim. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);SetNthreads(1);Alarm(300);tt:=Cputime();
R<T,x,y>:=PolynomialRing(GF(2),3);B:=FieldOfFractions(R);
F25:=eval Read("work/deg10/shared/mdmagma/models_X1_n/FFFc25.txt");
xx:=B!x;yy:=B!y;r:=(xx^2*yy-xx*yy+yy-1)/(xx*(xx*yy-1));s:=(xx*yy-yy+1)/(xx*yy);
b:=r*s*(r-1);c:=s*(r-1);
// The exact degree5 modular unit from the archived vanHoeij gonality table.
f3:=b;f4:=c;f5:=b-c;f6:=s-1;f7:=s-r;f8:=r*s-2*r+1;f9:=s^2-s-r+1;
f10:=xx-yy+1;f12:=xx-yy;f13:=xx^3*yy-xx^2*yy^2-xx^2*yy+xx*yy^2-yy+1;
h:=f3*f4*f6^3*f7^2*f10*f12^3/(f5*f8^2*f9^2*f13^3);
N:=R!Numerator(h);D:=R!Denominator(h);
print "N75_GON5_CANONICAL_RESULTANT_START",TotalDegree(N),TotalDegree(D),Cputime(tt);
res:=Resultant(F25,N-T*D,3);assert res ne 0;
K<t>:=RationalFunctionField(GF(2));S<X>:=PolynomialRing(K);proj:=hom<R->S|[S!t,X,S!0]>;
pol:=proj(res);factors:=Factorization(pol);
print "N75_GON5_CANONICAL_RESULTANT_FACTORS",[<Degree(z[1]),z[2]>:z in factors],Cputime(tt);
fs:=[z[1]:z in factors|Degree(z[1]) eq 5];assert #fs eq 1;
L<a>:=FunctionField(fs[1]);SY<Y>:=PolynomialRing(L);ev:=hom<R->SY|[SY!t,SY!a,Y]>;
gg:=GCD(ev(F25),ev(N-T*D));assert Degree(gg) eq 1;
ya:=-Coefficient(gg,0)/Coefficient(gg,1);assert Evaluate(ev(F25),ya) eq 0;
xa:=a;ra:=(xa^2*ya-xa*ya+ya-1)/(xa*(xa*ya-1));sa:=(xa*ya-ya+1)/(xa*ya);
ba:=ra*sa*(ra-1);ca:=sa*(ra-1);E:=EllipticCurve([1-ca,-ba,-ba,L!0,L!0]);P:=E![0,0];
assert 25*P eq E!0 and 5*P ne E!0 and Genus(L) eq 12;
print "N75_GON5_CANONICAL_BASE_MODEL_VERIFIED",Degree(L),Genus(L),Cputime(tt);
S3<W3>:=PolynomialRing(L);j25:=jInvariant(E);f:=W3^4+W3^3+1/j25;assert Degree(f) eq 4;
assert IsIrreducible(f);M<w>:=FunctionField(f);
print "N75_GON5_CANONICAL_DEGREE20_GENUS_START",Cputime(tt);
assert Genus(M) eq 73;
print "N75_GON5_CANONICAL_GENUS73_MODEL_VERIFIED",Cputime(tt);
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
Attach("work/deg10/shared/fast_hecke_11.m");
coords:=[M!t,M!a,w];
Attach("work/deg10/shared/fast_hecke.m");
load "work/degree9/cyclic/continuation_20261026/N75_geometry/noncyclic/deps/pair_af_place_linear_C48.m";
OL75:=MaximalOrderFinite(L);OM75:=MaximalOrderFinite(M);
n75,d75:=IntegralSplit(1/j25,OL75);assert M!d75*w in OM75;
function NativePoint75Relative(vals,wv)
 p0:=C48PairAFPlaceLinear([L!t,a],vals);d:=Degree(p0);F:=Parent(wv);assert d eq LCM([Degree(MinimalPolynomial(v,GF(2))):v in vals]);
 assert Valuation(j25,p0) eq 0 and &and[Valuation(c,p0) ge 0:c in Coefficients(f)|c ne 0] and Valuation(Discriminant(f),p0) eq 0;
 found:=false;
 for bits in [[1,0],[0,1],[1,1]] do theta:=&+[bits[i]*vals[i]:i in [1..2]];if Degree(MinimalPolynomial(theta,GF(2))) eq d then chosen:=bits;found:=true;break;end if;end for;assert found;
 e:=1;while wv^(2^(d*e)) ne wv do e+:=1;end while;assert e in [1,3];
 RF<V>:=PolynomialRing(F);hp:=&*[V-wv^(2^(d*k)):k in [0..e-1]];assert Degree(hp) eq e;
 mat:=Matrix(GF(2),[Eltseq(theta^i,GF(2)):i in [0..d-1]]);assert Rank(mat) eq d;thetaf:=chosen[1]*(L!t)+chosen[2]*a;hc:=[];
 for v in Coefficients(hp) do cs:=Eltseq(Solution(mat,Vector(GF(2),Eltseq(v,GF(2)))));hv:=&+[L!cs[i+1]*thetaf^i:i in [0..d-1]];assert Evaluate(hv,p0)^(2^d) eq Evaluate(hv,p0);Append(~hc,hv);end for;
 rel:=&+[M!hc[i+1]*w^i:i in [0..#hc-1]];ff,gg:=TwoGenerators(p0);assert ideal<Order(Ideal(p0))|ff,gg> eq Ideal(p0);
 IB:=ideal<OM75|OM75!(M!ff),OM75!(M!gg)>;I:=(IB+rel*OM75) meet ideal<OM75|1>;
 printf "DECIC75_RELATIVE_POINT_IDEAL base_degree=%o residue_degree=%o time=%o\n",d,e,Cputime(tt);
 assert IsIntegral(I) and IsPrime(I);pp:=Place(I);assert Degree(pp) eq d*e;
 assert Valuation(M!ff,pp) gt 0 and Valuation(M!gg,pp) gt 0 and Valuation(rel,pp) gt 0 and Valuation(M!UniformizingElement(p0),pp) eq 1;
 assert Degree(pp) eq LCM([Degree(MinimalPolynomial(v,GF(2))):v in vals cat [wv]]);
 return pp;
end function;

function AbsoluteRecord(vals)
 F:=Universe(vals);d:=LCM([Degree(MinimalPolynomial(v,GF(2))):v in vals]);assert d in [10,30];
 found:=false;
 for v in vals do if Degree(MinimalPolynomial(v,GF(2))) eq d then theta:=v;found:=true;break;end if;end for;
 assert found;mp:=MinimalPolynomial(theta,GF(2));
 mat:=Matrix(GF(2),[Eltseq(theta^j,GF(2)):j in [0..d-1]]);assert Rank(mat) eq d;
 enc:=[];
 for v in vals do cc:=Eltseq(Solution(mat,Vector(GF(2),Eltseq(v,GF(2)))));Append(~enc,&+[Integers()!cc[j+1]*2^j:j in [0..d-1]]);end for;
 return [Integers()!c:c in Coefficients(mp)],enc;
end function;
rowsfile:=Sprintf("work/deg10/continuation_20261030b/C75/moduli_job%o.txt",JobIndex);
fo:=Open(rowsfile,"w");delete fo;
function ModuliY(E,P25,P3,K,coords)
 assert Order(P25) eq 25 and Order(P3) eq 3;
 a1,a2,a3,a4,a6:=Explode(aInvariants(E));x0:=P25[1];y0:=P25[2];
 aa3:=2*y0+a3+a1*x0;aa4:=3*x0^2+2*x0*a2+a4-a1*y0;sl:=aa4/aa3;
 aa2:=3*x0+a2-sl*a1-sl^2;sc:=aa2/aa3;
 b:=-aa2^3/aa3^2;c:=1-(a1+2*sl)*sc;
 EC:=EllipticCurve([1-c,-b,-b,BaseRing(E)!0,BaseRing(E)!0]);
 RR:=EC![sc^2*(P3[1]-x0),sc^3*(P3[2]-y0-sl*(P3[1]-x0))];assert Order(RR) eq 3;
 aa:=aInvariants(EC);a1,a2,a3,a4,a6:=Explode(aa);
 xc:=(RR[1]-a3/a1)/a1^2;assert xc^4+xc^3+1/jInvariant(E) eq 0;
 xy:=TateNormalForm_xy(E,P25);xx,yy:=Explode(xy);
 r:=(xx^2*yy-xx*yy+yy-1)/(xx*(xx*yy-1));s:=(xx*yy-yy+1)/(xx*yy);
 f3:=b;f4:=c;f5:=b-c;f6:=s-1;f7:=s-r;f8:=r*s-2*r+1;f9:=s^2-s-r+1;
 f10:=xx-yy+1;f12:=xx-yy;f13:=xx^3*yy-xx^2*yy^2-xx^2*yy+xx*yy^2-yy+1;
 hv:=f3*f4*f6^3*f7^2*f10*f12^3/(f5*f8^2*f9^2*f13^3);
 pp:=NativePoint75Relative([hv,xx],xc);
 poly,enc:=AbsoluteRecord(aInvariants(E) cat [P25[1],P25[2],P3[1],P3[2],b,c,xx,yy,hv,xc]);
 row:=[*poly,enc,Degree(pp)*];
 Write(rowsfile,Sprint(row));return pp;
end function;

function ER(v)
 return [[Integers()!c:c in Coefficients(Numerator(v))],[Integers()!c:c in Coefficients(Denominator(v))]];
end function;
function EL(v) return [ER(c):c in Eltseq(L!v)];end function;
function EM(v) return [EL(c):c in Eltseq(M!v)];end function;
function PK(p)
 J:=Ideal(p);ff,gg:=TwoGenerators(p);assert ideal<Order(J)|ff,gg> eq J;return [*IsFinite(p) select 1 else 0,Degree(p),[EM(ff),EM(gg)]*];
end function;
function PD(D)
 ps,es:=Support(D);return [*[*PK(ps[k]),es[k]*]:k in [1..#ps]*];
end function;
load "work/deg10/continuation_20261030b/C75/census10_records.m";
FF<aa>:=ext<GF(2)|PolynomialRing(GF(2))!FieldPolynomial>;evp:=func<v|Evaluate(PolynomialRing(GF(2))!v,aa)>;
jpolys:=Sort(Setseq({rr[6]:rr in Records}));assert #jpolys eq 6;
for ci in CaseIndices do
 rows:=[rr:rr in Records|rr[6] eq jpolys[ci]];assert #rows eq 10;rr:=rows[1];EC:=EllipticCurve([evp(v):v in rr[2]]);AG,am:=AbelianGroup(EC);assert #EC eq 1050 and Invariants(AG) eq [1050];P75:=14*am(AG.1);assert Order(P75) eq 75;P25:=3*P75;P3:=25*P75;
 Psource:=ModuliY(EC,P25,P3,M,coords);assert Degree(Psource) eq 10;
 printf "DECIC75_INPUT_READY case=%o degree10=true time=%o\n",ci,Cputime(tt);
 q:=7;terms:=MDIsogeniesFastPoints(EC,[P25,P3],q);H:=DivisorGroup(M)!0;fastmeta:=[];
 for term in terms do
  EG,points,e:=Explode(term);pp:=ModuliY(EG,points[1],points[2],M,coords);mult:=e*10/Degree(pp);assert Denominator(mult) eq 1;H+:=(Integers()!mult)*Divisor(pp);Append(~fastmeta,[e,Degree(pp),Integers()!mult]);
 end for;
 assert Degree(H) eq 80;delta:=ModuliY(EC,q*P25,q*P3,M,coords);assert Degree(delta) eq 10;Dtest:=H-Divisor(Psource)-7*Divisor(delta);assert Degree(Dtest) eq 0;
 packet:=[*ci,[ER(c):c in Coefficients(DefiningPolynomial(L))],[EL(c):c in Coefficients(DefiningPolynomial(M))],EL(j25),PD(Divisor(Psource)),PD(Divisor(delta)),PD(H),PD(Dtest),fastmeta,rr*];
 Write(Sprintf("work/deg10/continuation_20261030b/C75/proof_packet_case%o.m",ci),Sprint(packet):Overwrite:=true);
 bb:=Basis(Dtest:Reduction:=false,Simplification:="None");assert #bb eq 0;
 printf "DECIC75_A7_RR_COMPLETE case=%o h0=0 Tdegree80=true point_convention=true metadata=%o time=%o\n",ci,fastmeta,Cputime(tt);
 native:=MDIsogenies(EC,7);assert #native eq 8;nativeplaces:=[];
 for phi in native do
  EE:=Domain(phi);Q25:=phi(EE![BaseRing(EE)!v:v in Eltseq(P25)]);Q3:=phi(EE![BaseRing(EE)!v:v in Eltseq(P3)]);qq:=ModuliY(Codomain(phi),Q25,Q3,M,coords);Append(~nativeplaces,qq);
 end for;
 common:=LCM([Degree(qq):qq in nativeplaces]);assert common mod 10 eq 0;ND:=&+[DivisorGroup(M)|(common div Degree(qq))*Divisor(qq):qq in nativeplaces];assert ND eq (common div 10)*H;
 printf "DECIC75_NATIVE_GEOMETRIC_HECKE_CONTROL_PASS case=%o common=%o geometric_kernels=8 time=%o\n",ci,common,Cputime(tt);
end for;
printf "DECIC75_NC10_JOB_COMPLETE job=%o cases=%o time=%o\n",JobIndex,CaseIndices,Cputime(tt);quit;
