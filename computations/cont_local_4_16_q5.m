SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(5*10^9);t0:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
Attach("work/deg10/continuation_20261029/deps/fast_hecke.m");
Attach("work/deg10/continuation_20261029/deps/fast_hecke_11.m");
X:=MDX11(4,16,GF(5));
C:=Curve(X);assert Genus(C) eq 13;
cs:=Cusps(X);ds:=[Degree(c):c in cs];
print "DECIC_LOCAL_CUSPS",4,16,5,ds,Cputime(t0);
ns:=[];
for d in [1..5] do xs:=NoncuspidalPlaces(X,d);Append(~ns,#xs);print "DECIC_LOCAL_COUNT",4,16,5,d,#xs,Cputime(t0);end for;
J:=JH(64,[17]);assert Dimension(J) eq 13;
p:=Characteristic(GF(5));e:=Degree(GF(5));f:=FrobeniusPolynomial(J,p);assert Degree(f) eq 26;tr:=[];
for i in [1..5*e] do
 a:=-i*Coefficient(f,26-i);
 for j in [1..i-1] do a-:=Coefficient(f,26-j)*tr[i-j];end for;
 Append(~tr,Integers()!a);
end for;
cl:=[];
for d in [1..5] do
 a:=5^d+1-tr[e*d];for j in Divisors(d) do if j lt d then a-:=j*cl[j];end if;end for;assert a mod d eq 0;Append(~cl,a div d);
end for;
assert cl eq [ns[d]+#[x:x in ds|x eq d]:d in [1..5]];
Write("work/deg10/continuation_20261029/cont_local_4_16_q5_census.m",Sprint([*4,16,5,ds,ns,Coefficients(f),cl*]));
print "DECIC_LOCAL_CENSUS_COMPLETE",4,16,5,ns,Cputime(t0);
AF:=AlgorithmicFunctionField(FunctionField(C));c0:=[c:c in cs|Degree(c) eq 1][1];
G,_,phi:=ClassGroup(AF:ReductionDivisor:=FunctionFieldDivisor(Divisor(c0)),Proof:=true);
print "DECIC_LOCAL_CLASSGROUP_PROVEN",Invariants(G),Cputime(t0);
quit;
