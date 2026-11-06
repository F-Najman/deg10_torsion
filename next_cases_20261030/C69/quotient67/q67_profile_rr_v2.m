SetColumns(0);SetNthreads(1);SetQuitOnError(true);SetMemoryLimit(16*10^9);
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");Attach("work/degree9/shared/fast_hecke.m");
Alarm(300);tt:=Cputime();record:=eval Read("work/deg10/next_cases_20261030/C69/direct_full/selected_input.m");assert record[1..3] eq [*69,10,2*];
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
// The order-three Tate family identifies the quartic with X0(3) over j.
F0<t>:=FunctionField(GF(2));E3:=EllipticCurve([F0!1,F0!0,t,F0!0,F0!0]);P3:=E3![0,0];assert Order(P3) eq 3 and jInvariant(E3) eq 1/(t^3*(t+1));
RV<v>:=PolynomialRing(K);f:=v^4+v^3+1/jK;assert IsIrreducible(f);Z<s>:=FunctionField(f);assert Genus(Z) eq 67;
printf "Q67_QUARTIC_MODEL_READY genus=67 irreducible=true time=%o\n",Cputime(tt);
// The simple root v=1 at a rational cusp proves the exact constant field.
rc:=[P:P in Poles(jK)|Degree(P) eq 1];assert #rc eq 11;
assert Valuation(Evaluate(f,K!1),rc[1]) gt 0 and Valuation(Evaluate(Derivative(f),K!1),rc[1]) eq 0;
printf "Q67_V2_CONSTANTS_HENSEL_VERIFIED rational_simple_cusp_root=true time=%o\n",Cputime(tt);
out:=Open("work/deg10/next_cases_20261030/C69/quotient67/model_profile_v2.json","w");fprintf out,"{\"base_genus\":12,\"quotient_genus\":67,\"relative_degree\":4,\"irreducible\":true,\"constants_F2_by_simple_cusp_root\":true,\"model_profile_complete\":true,\"cpu_seconds\":%o}\n",Cputime(tt);delete out;
Alarm(1200);rt:=Cputime();
NumJ,DenJ:=IntegralSplit(jK,O0);NumJ:=K!NumJ;DenJ:=K!DenJ;
OZ:=MaximalOrderFinite(Z);IZ:=(Z!NumJ)*s;assert IZ^4+(Z!NumJ)*IZ^3+(Z!(NumJ^3*DenJ)) eq 0;assert OZ!IZ in OZ;
printf "Q67_V2_FINITE_ORDER_READY time=%o\n",Cputime(tt);
function UniqueNCLift(P,gen)
 d:=Degree(P);assert Valuation(NumJ,P) eq 0 and Valuation(DenJ,P) eq 0;
 L:=ResidueClassField(P);RL<w>:=PolynomialRing(L);ff:=RL![Evaluate(c,P):c in Coefficients(f)];fac:=Factorization(ff);
 assert Sort([Degree(c[1]):c in fac]) eq [1,3] and &and[c[2] eq 1:c in fac];roots:=Roots(ff);assert #roots eq 1;beta:=roots[1][1];
 gval:=Evaluate(gen,P);assert Degree(MinimalPolynomial(gval,GF(2))) eq d;MB:=Matrix(GF(2),[Eltseq(gval^i,GF(2)):i in [0..d-1]]);assert Rank(MB) eq d;
 co:=Solution(MB,Vector(GF(2),Eltseq(beta,GF(2))));lift:=&+[K!co[k+1]*gen^k:k in [0..d-1]];assert Evaluate(lift,P) eq beta and O0!lift in O0;
 gg,hh:=TwoGenerators(P);assert ideal<Order(Ideal(P))|gg,hh> eq Ideal(P);u:=UniformizingElement(P);
 J:=ideal<OZ|OZ!(Z!gg),OZ!(Z!hh),OZ!(IZ-(Z!NumJ)*(Z!lift))>;assert IsPrime(J);Q:=Place(J);assert Degree(Q) eq d;
 assert Valuation(Z!gg,Q) gt 0 and Valuation(Z!hh,Q) gt 0 and Valuation(Z!u,Q) eq 1 and Valuation(s-Z!lift,Q) gt 0;
 return Q;
end function;
BP,genBP:=BasePlace(E,R23);assert Degree(BP) eq 10;PZ:=UniqueNCLift(BP,genBP);assert Degree(PZ) eq 10;
printf "Q67_V2_ACTUAL_INPUT_LIFT_READY base_case=2 degree10=true unique_linear_factor=true time=%o\n",Cputime(tt);
T:=DivisorGroup(Z)!0;imgs:=[* *];baseplaces:=[* *];coverplaces:=[* *];isogs:=MDIsogeniesFast(E,PN,5);assert &+[q[3]:q in isogs] eq 6;
for idx in [1..#isogs] do
 EE,PP,w:=Explode(isogs[idx]);assert Order(PP) eq 69;Q23:=3*PP;assert Order(Q23) eq 23;B,genB:=BasePlace(EE,Q23);d:=Degree(B);
 xy69:=TateNormalForm_xy(EE,PP);d69:=LCM([Degree(MinimalPolynomial(c,GF(2))):c in xy69]);dj:=Degree(MinimalPolynomial(jInvariant(EE),GF(2)));assert d eq d69 and dj eq d and d in [10,50];
 QZ:=UniqueNCLift(B,genB);assert Degree(QZ) eq d and (10*w) mod d eq 0;mult:=(10*w) div d;T+:=mult*Divisor(QZ);Append(~imgs,[*idx,d,d69,dj,w,mult*]);Append(~baseplaces,B);Append(~coverplaces,QZ);
 printf "Q67_HEKCE_ORBIT_RESTRICTED index=%o degree=%o full69degree=%o jdegree=%o geometric_observations=%o norm_weight=%o time=%o\n",idx,d,d69,dj,w,mult,Cputime(tt);
end for;
assert Degree(T) eq 60;BD,genBD:=BasePlace(E,5*R23);assert Degree(BD) eq 10;DP:=UniqueNCLift(BD,genBD);A5:=T-Divisor(PZ)-5*Divisor(DP);assert Degree(A5) eq 0;
printf "Q67_A5_DIVISOR_READY degree=0 hecke_degree=60 diamond5_retained=true time=%o\n",Cputime(tt);
function ER(v)
 return [[Integers()!c:c in Coefficients(Numerator(v))],[Integers()!c:c in Coefficients(Denominator(v))]];
end function;
function EK(v) return [ER(c):c in Eltseq(K!v)];end function;
function EZ(v) return [EK(c):c in Eltseq(Z!v)];end function;
function DR(v)
 RR:=PolynomialRing(GF(2));return (BaseField(K)!(RR!v[1]))/(BaseField(K)!(RR!v[2]));
end function;
function DK(v) return &+[K!DR(v[i+1])*K.1^i:i in [0..#v-1]];end function;
function DZ(v) return &+[Z!DK(v[i+1])*Z.1^i:i in [0..#v-1]];end function;
function PK(P)
 J:=Ideal(P);g,h:=TwoGenerators(P);assert ideal<Order(J)|g,h> eq J;v:=[EK(g),EK(h)];assert DK(v[1]) eq g and DK(v[2]) eq h;
 return [*IsFinite(P) select 1 else 0,Degree(P),v*];
end function;
function PZpack(P)
 J:=Ideal(P);g,h:=TwoGenerators(P);assert ideal<Order(J)|g,h> eq J;v:=[EZ(g),EZ(h)];assert DZ(v[1]) eq g and DZ(v[2]) eq h;
 return [*IsFinite(P) select 1 else 0,Degree(P),v*];
end function;
function PD(D)
 ps,es:=Support(D);return [*[*PZpack(ps[i]),es[i]*]:i in [1..#ps]*];
end function;
function DD(rows)
 D:=DivisorGroup(Z)!0;
 for row in rows do p,e:=Explode(row);OO:=p[1] eq 1 select MaximalOrderFinite(Z) else MaximalOrderInfinite(Z);J:=ideal<OO|[DZ(v):v in p[3]]>;assert IsPrime(J);Q:=Place(J);assert Degree(Q) eq p[2];D+:=e*Divisor(Q);end for;return D;
end function;
kp:=[ER(c):c in Coefficients(DefiningPolynomial(K))];zp:=[EK(c):c in Coefficients(DefiningPolynomial(Z))];assert [DR(c):c in kp] eq Coefficients(DefiningPolynomial(K)) and [DK(c):c in zp] eq Coefficients(DefiningPolynomial(Z));
Dsource:=PD(Divisor(PZ));Ddiamond:=PD(Divisor(DP));Dhecke:=PD(T);Dtest:=PD(A5);
assert DD(Dsource) eq Divisor(PZ) and DD(Ddiamond) eq Divisor(DP) and DD(Dhecke) eq T and DD(Dtest) eq A5;
pairs:=[*[*PK(BP),PZpack(PZ),10,1,1*], [*PK(BD),PZpack(DP),10,1,1*]*];
for i in [1..#baseplaces] do Append(~pairs,[*PK(baseplaces[i]),PZpack(coverplaces[i]),imgs[i][2],imgs[i][5],imgs[i][6]*]);end for;
packet:=[*"Q67_PACKET_V2",kp,zp,EK(jK),Dsource,Ddiamond,Dhecke,Dtest,imgs,pairs,record*];
Write("work/deg10/next_cases_20261030/C69/quotient67/proof_packet_v2.m",Sprint(packet):Overwrite:=true);
printf "Q67_V2_PROOF_PACKET_COMPLETE two_generators=true ideal_equalities=true function_roundtrips=true divisor_roundtrips=true time=%o\n",Cputime(tt);
principal:=IsPrincipal(A5);out:=Open("work/deg10/next_cases_20261030/C69/quotient67/rr_result_v2.json","w");fprintf out,"{\"N\":69,\"quotient_genus\":67,\"base_case\":2,\"point_degree\":10,\"operator\":\"T5-1-5<5>\",\"principal\":%o,\"actual_source_lift_and_all_Hecke_lifts_unique\":true,\"isogeny_orbits\":%o,\"completed\":true,\"rr_cpu_seconds\":%o,\"total_cpu_seconds\":%o}\n",principal select "true" else "false",[[Integers()!v:v in r]:r in imgs],Cputime(rt),Cputime(tt);delete out;
printf "Q67_V2_A5_RR_COMPLETE base_case=2 principal=%o unique_actual_lifts=true rr_time=%o total_time=%o\n",principal,Cputime(rt),Cputime(tt);quit;
