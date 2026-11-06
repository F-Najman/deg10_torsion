SetColumns(0);SetNthreads(1);SetQuitOnError(true);SetMemoryLimit(16*10^9);
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");Attach("work/degree9/shared/fast_hecke.m");
Alarm(600);tt:=Cputime();record:=eval Read("work/deg10/next_cases_20261030/C69/direct_full/selected_input.m");assert record[1..3] eq [*69,10,2*];
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
 return P;
end function;
// The order-three Tate family identifies the quartic with X0(3) over j.
F0<t>:=FunctionField(GF(2));E3:=EllipticCurve([F0!1,F0!0,t,F0!0,F0!0]);P3:=E3![0,0];assert Order(P3) eq 3 and jInvariant(E3) eq 1/(t^3*(t+1));
RV<v>:=PolynomialRing(K);f:=v^4+v^3+1/jK;assert IsIrreducible(f);Z<s>:=FunctionField(f);assert Genus(Z) eq 67;
printf "Q67_QUARTIC_MODEL_READY genus=67 irreducible=true time=%o\n",Cputime(tt);
function RestrictedPlaces(P)
 assert BaseField(Z) cmpeq FunctionField(P);g,h:=TwoGenerators(P);u:=UniformizingElement(P);assert Valuation(g,P) gt 0 and Valuation(h,P) gt 0 and Valuation(u,P) eq 1;
 raw:=Decomposition(Z,P);good:=[Q:Q in raw|Valuation(Z!g,Q) gt 0 and Valuation(Z!h,Q) gt 0];assert #good gt 0;
 assert &+[Degree(Q)*Valuation(Z!u,Q):Q in good] eq 4*Degree(P);return good,[Valuation(Z!u,Q):Q in good];
end function;
// A degree-one cusp verifies the exact constant field of this component.
rc:=Places(K,1);assert #rc eq 11 and &and[Valuation(jK,P) lt 0:P in rc];dc,es:=RestrictedPlaces(rc[1]);assert #dc eq 2 and &and[Degree(Q) eq 1:Q in dc] and Sort(es) eq [1,3];
printf "Q67_RATIONAL_CUSP_CONSTANTS_VERIFIED two_rational_lifts=true ramification=1,3 time=%o\n",Cputime(tt);
function UniqueNCLift(P)
 dd,es:=RestrictedPlaces(P);d:=Degree(P);assert #dd eq 2 and Sort([Degree(Q):Q in dd]) eq [d,3*d] and es eq [1,1];same:=[Q:Q in dd|Degree(Q) eq d];assert #same eq 1;return same[1];
end function;
BP:=BasePlace(E,R23);assert Degree(BP) eq 10;PZ:=UniqueNCLift(BP);assert Degree(PZ) eq 10;
out:=Open("work/deg10/next_cases_20261030/C69/quotient67/profile.json","w");fprintf out,"{\"base_genus\":12,\"quotient_genus\":67,\"relative_degree\":4,\"base_case\":2,\"input_degree\":10,\"rational_cusp_constants_verified\":true,\"actual_input_unique_lift\":true,\"profile_complete\":true,\"cpu_seconds\":%o}\n",Cputime(tt);delete out;
printf "Q67_PROFILE_COMPLETE base_case=2 input_degree=10 unique_lift=true genus=67 time=%o\n",Cputime(tt);
Alarm(1200);rt:=Cputime();T:=DivisorGroup(Z)!0;imgs:=[* *];isogs:=MDIsogeniesFast(E,PN,5);assert &+[q[3]:q in isogs] eq 6;
for idx in [1..#isogs] do
 EE,PP,w:=Explode(isogs[idx]);assert Order(PP) eq 69;Q23:=3*PP;assert Order(Q23) eq 23;B:=BasePlace(EE,Q23);d:=Degree(B);
 xy69:=TateNormalForm_xy(EE,PP);d69:=LCM([Degree(MinimalPolynomial(c,GF(2))):c in xy69]);dj:=Degree(MinimalPolynomial(jInvariant(EE),GF(2)));assert d eq d69 and dj eq d and d in [10,50];
 QZ:=UniqueNCLift(B);assert Degree(QZ) eq d and (10*w) mod d eq 0;mult:=(10*w) div d;T+:=mult*Divisor(QZ);Append(~imgs,[*idx,d,d69,dj,w,mult*]);
 printf "Q67_HEKCE_ORBIT_RESTRICTED index=%o degree=%o full69degree=%o jdegree=%o geometric_observations=%o norm_weight=%o time=%o\n",idx,d,d69,dj,w,mult,Cputime(tt);
end for;
assert Degree(T) eq 60;BD:=BasePlace(E,5*R23);assert Degree(BD) eq 10;DP:=UniqueNCLift(BD);A5:=T-Divisor(PZ)-5*Divisor(DP);assert Degree(A5) eq 0;
printf "Q67_A5_DIVISOR_READY degree=0 hecke_degree=60 diamond5_retained=true time=%o\n",Cputime(tt);
principal:=IsPrincipal(A5);out:=Open("work/deg10/next_cases_20261030/C69/quotient67/rr_result.json","w");fprintf out,"{\"N\":69,\"quotient_genus\":67,\"base_case\":2,\"point_degree\":10,\"operator\":\"T5-1-5<5>\",\"principal\":%o,\"actual_source_lift_and_all_Hecke_lifts_unique\":true,\"isogeny_orbits\":%o,\"completed\":true,\"rr_cpu_seconds\":%o,\"total_cpu_seconds\":%o}\n",principal select "true" else "false",[[Integers()!v:v in r]:r in imgs],Cputime(rt),Cputime(tt);delete out;
printf "Q67_A5_RR_COMPLETE base_case=2 principal=%o unique_actual_lifts=true rr_time=%o total_time=%o\n",principal,Cputime(rt),Cputime(tt);quit;
