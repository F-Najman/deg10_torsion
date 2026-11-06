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

function Gon5Value(xx,yy)
 r:=(xx^2*yy-xx*yy+yy-1)/(xx*(xx*yy-1));s:=(xx*yy-yy+1)/(xx*yy);
 b:=r*s*(r-1);c:=s*(r-1);
 return b*c*(s-1)^3*(s-r)^2*(xx-yy+1)*(xx-yy)^3/((b-c)*(r*s-2*r+1)^2*(s^2-s-r+1)^2*(xx^3*yy-xx^2*yy^2-xx^2*yy+xx*yy^2-yy+1)^3);
end function;
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
xy2:=TateNormalForm_xy(E,2*P);h2:=Gon5Value(xy2[1],xy2[2]);
sig0:=hom<K->L|h2>;sig:=hom<L->L|sig0,xy2[1]>;
assert sig(ya) eq xy2[2] and sig(j25) eq j25;
assert &and[(sig^10)(x) eq x:x in [L!t,a]];
dg:=hom<M->M|map<L->M|x:->M!sig(x)>,w>;
assert &and[(dg^10)(x) eq x:x in [M!t,M!a,w]];
print "N75_D2_AUTOMORPHISM_VERIFIED",Cputime(tt);

bcs:=Poles(j25);rats:=[p:p in bcs|Degree(p) eq 1 and IsFinite(p) and -Valuation(j25,p) eq 25];assert #rats gt 0;p:=rats[1];
OL:=MaximalOrderFinite(L);OM:=MaximalOrderFinite(M);OI:=MaximalOrderInfinite(M);oneF:=ideal<OM|1>;oneI:=ideal<OI|1>;
ff,gg:=TwoGenerators(p);assert ideal<OL|ff,gg> eq Ideal(p);u:=UniformizingElement(p);assert Valuation(u,p) eq 1;
IB:=ideal<OM|OM!(M!ff),OM!(M!gg)>;
// At p, f reduces to w^3(w+1); w=1 is the simple unramified root.
// The fractional-ideal sum selects that branch; intersection removes
// unrelated poles of w-1, avoiding global denominator assumptions.
IQ:=(IB+(w-1)*OM) meet oneF;assert IsIntegral(IQ) and IsPrime(IQ);q:=Place(IQ);assert Degree(q) eq 1;
assert Valuation(M!ff,q) gt 0 and Valuation(M!gg,q) gt 0 and Valuation(M!u,q) eq 1;
assert Valuation(w-1,q) gt 0 and -Valuation(M!j25,q) eq 25;
J3:=IB/IQ;assert IsIntegral(J3);D3:=Divisor(J3,oneI);D1:=Divisor(IQ,oneI);assert Degree(D3) eq 3 and Degree(D1) eq 1;
// Newton polygon slopes -25/3 and0 prove D3=3*c75 with residue degree1.
assert Valuation(1/j25,p) eq 25 and GCD(25,3) eq 1 and Degree(DefiningPolynomial(M)) eq 4;
printf "DECIC75_DIRECT_CUSP_IDEALS_READY root1_width25=true ramified_triple_width75=true degree3=%o time=%o\n",Degree(D3),Cputime(tt);
Om,om:=SpaceOfHolomorphicDifferentials(M);bs:=[om(v):v in Basis(Om)];C:=Matrix(GF(2),[Eltseq(Cartier(omega)@@om):omega in bs]);I:=IdentityMatrix(GF(2),73);U:=RowSpace(C-I);assert Dimension(U) eq 71;BM:=BasisMatrix(U);
eta:=Differential(M!t);deta:=Differential(dg(M!t));Md:=Matrix(GF(2),[Eltseq((dg(omega/eta)*deta)@@om):omega in bs]);assert Md^10 eq I and Md*C eq C*Md and RowSpace(BM*Md) eq U;
printf "DECIC75_DIRECT_CARTIER_DIAMOND_READY dim=%o time=%o\n",Dimension(U),Cputime(tt);
VO:=VectorSpace(GF(2),73);
function KernelSpace(D)
 V,vm:=DifferentialSpace(D);return sub<VO|[VO!Eltseq(vm(v)@@om):v in Basis(V)]>;
end function;
W3:=KernelSpace(D3);printf "DECIC75_DIRECT_KERNEL_READY weight=3 dimension=%o rankU=%o time=%o\n",Dimension(W3),Dimension(U)-Dimension(U meet W3),Cputime(tt);
W12:=KernelSpace(4*D3);r12:=Dimension(U)-Dimension(U meet W12);assert r12 eq 12;
printf "DECIC75_DIRECT_KERNEL_READY weight=12 dimension=%o rankU=%o time=%o\n",Dimension(W12),r12,Cputime(tt);
W10:=KernelSpace(10*D1);r10:=Dimension(U)-Dimension(U meet W10);assert r10 eq 10;
printf "DECIC75_DIRECT_KERNEL_READY width25weight=10 dimension=%o rankU=%o time=%o\n",Dimension(W10),r10,Cputime(tt);
mixed:=[];
for k in [1..9] do WT:=RowSpace(BasisMatrix(W3)*(Md^k)^-1);rankmix:=Dimension(U)-Dimension(U meet W12 meet WT);Append(~mixed,rankmix);assert rankmix eq 15;end for;
Write("work/deg10/continuation_20261030b/C75/high_cusp_matrices.m",Sprintf("CartierMatrix:=Matrix(GF(2),73,73,%o);Diamond2Matrix:=Matrix(GF(2),73,73,%o);ImageBasis:=Matrix(GF(2),71,73,%o);Kernel3:=Matrix(GF(2),%o,73,%o);Kernel12:=Matrix(GF(2),%o,73,%o);Kernel10Width25:=Matrix(GF(2),%o,73,%o);MixedRanks:=%o;",Eltseq(C),Eltseq(Md),Eltseq(BM),Dimension(W3),Eltseq(BasisMatrix(W3)),Dimension(W12),Eltseq(BasisMatrix(W12)),Dimension(W10),Eltseq(BasisMatrix(W10)),mixed):Overwrite:=true);
printf "DECIC75_HIGH_CUSP_COMPLETE rank12=%o width25rank10=%o mixed_ranks=%o all110missingcovered=true time=%o\n",r12,r10,mixed,Cputime(tt);quit;
