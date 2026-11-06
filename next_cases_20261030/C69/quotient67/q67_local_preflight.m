SetColumns(0);SetNthreads(1);SetQuitOnError(true);SetMemoryLimit(16*10^9);
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");Attach("work/degree9/shared/fast_hecke.m");
Alarm(120);tt:=Cputime();record:=eval Read("work/deg10/next_cases_20261030/C69/direct_full/selected_input.m");assert record[1..3] eq [*69,10,2*];
RP<z>:=PolynomialRing(GF(2));F<a>:=ext<GF(2)|RP!record[5]>;fromvec:=func<v|&+[F!v[k]*a^(k-1):k in [1..#v]]>;
E:=EllipticCurve([fromvec(v):v in record[6]]);PN:=E![fromvec(v):v in record[7]];assert Order(PN) eq 69 and #E eq 966;R23:=3*PN;assert Order(R23) eq 23;
X23:=MDX1(23,GF(2));C23:=Curve(X23);K,am:=AlgorithmicFunctionField(FunctionField(C23));assert Genus(K) eq 12;
jK:=am(jInvariant(EllipticCurve(X23`_E)));assert Degree(jK) eq 264;
O0:=MaximalOrderFinite(K);IX:=K!BaseField(K).1;numY,denY:=IntegralSplit(K.1,O0);IY:=K!numY;denK:=K!denY;assert &and[c eq 0:c in Eltseq(denK)[2..Degree(K)]];
fx:=IX@@am;fy:=IY@@am;fd:=denK@@am;assert O0!IX in O0 and O0!IY in O0;
printf "Q67_BASE_X123_READY genus=12 jdegree=264 time=%o\n",Cputime(tt);
function BasePlace(E,Q)
 assert Order(Q) eq 23;L:=BaseRing(E);xy:=TateNormalForm_xy(E,Q);pt:=C23(L)!xy;assert Evaluate(fd,pt) ne 0;
 vals:=[Evaluate(fx,pt),Evaluate(fy,pt)];degrees:=[Degree(MinimalPolynomial(v,GF(2))):v in vals];d:=LCM(degrees);
 candidates:=[<IX,vals[1]>,<IY,vals[2]>,<IX+IY,vals[1]+vals[2]>];found:=false;
 for q in candidates do mp:=MinimalPolynomial(q[2],GF(2));if Degree(mp) eq d then gen:=q[1];value:=q[2];found:=true;break;end if;end for;assert found;
 MB:=Matrix(GF(2),[Eltseq(value^i,GF(2)):i in [0..d-1]]);assert Rank(MB) eq d;gens:=[O0!Evaluate(mp,gen)];
 for i in [1,2] do co:=Solution(MB,Vector(GF(2),Eltseq(vals[i],GF(2))));pol:=&+[K!co[k+1]*gen^k:k in [0..d-1]];Append(~gens,O0!([IX,IY][i]-pol));end for;
 I:=ideal<O0|gens>;assert IsPrime(I);P:=Place(I);assert Degree(P) eq d;
 assert Valuation(jK,P) ge 0 and MinimalPolynomial(Evaluate(jK,P),GF(2)) eq MinimalPolynomial(jInvariant(E),GF(2));
 return P,gen;
end function;

NumJ,DenJ:=IntegralSplit(jK,O0);NumJ:=K!NumJ;DenJ:=K!DenJ;
RV<v>:=PolynomialRing(K);f:=v^4+v^3+1/jK;
rc:=[P:P in Poles(jK)|Degree(P) eq 1];assert #rc eq 11 and Valuation(Evaluate(f,K!1),rc[1]) gt 0 and Valuation(Evaluate(Derivative(f),K!1),rc[1]) eq 0;
printf "Q67_PREFLIGHT_SIMPLE_CUSP_ROOT_PASS time=%o\n",Cputime(tt);
function TestFactor(E,Q,label)
 BP,gen:=BasePlace(E,Q);d:=Degree(BP);assert Valuation(NumJ,BP) eq 0 and Valuation(DenJ,BP) eq 0;
 L:=ResidueClassField(BP);RPV<w>:=PolynomialRing(L);ff:=RPV![Evaluate(c,BP):c in Coefficients(f)];fac:=Factorization(ff);
 assert Sort([Degree(c[1]):c in fac]) eq [1,3] and &and[c[2] eq 1:c in fac];
 printf "Q67_PREFLIGHT_UNIQUE_LINEAR_FACTOR label=%o degree=%o factors=1,3 integral_multiplier_unit=true time=%o\n",label,d,Cputime(tt);
 return d;
end function;
d:=TestFactor(E,R23,"input");assert d eq 10;
isogs:=MDIsogeniesFast(E,PN,5);for i in [1..#isogs] do EE,PP,w:=Explode(isogs[i]);d:=TestFactor(EE,3*PP,Sprintf("hecke%o",i));assert d eq 10*w and d in [10,50];end for;
d:=TestFactor(E,5*R23,"diamond5");assert d eq 10;
print "Q67_LOCAL_PREFLIGHT_COMPLETE";quit;
