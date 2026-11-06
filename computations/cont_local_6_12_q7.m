SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(5*10^9);t0:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
Attach("work/deg10/continuation_20261029/deps/fast_hecke.m");
Attach("work/deg10/continuation_20261029/deps/fast_hecke_11.m");
Attach("work/deg10/continuation_20261029/deps/model_6_12.m");
X:=MDX612(GF(7));
C:=Curve(X);assert Genus(C) eq 9;
cs:=Cusps(X);ds:=[Degree(c):c in cs];
print "DECIC_LOCAL_CUSPS",6,12,7,ds,Cputime(t0);
ns:=[];
for d in [1..5] do xs:=NoncuspidalPlaces(X,d);Append(~ns,#xs);print "DECIC_LOCAL_COUNT",6,12,7,d,#xs,Cputime(t0);end for;
J:=JH(72,[13]);assert Dimension(J) eq 9;
p:=Characteristic(GF(7));e:=Degree(GF(7));f:=FrobeniusPolynomial(J,p);assert Degree(f) eq 18;tr:=[];
for i in [1..5*e] do
 a:=-i*Coefficient(f,18-i);
 for j in [1..i-1] do a-:=Coefficient(f,18-j)*tr[i-j];end for;
 Append(~tr,Integers()!a);
end for;
cl:=[];
for d in [1..5] do
 a:=7^d+1-tr[e*d];for j in Divisors(d) do if j lt d then a-:=j*cl[j];end if;end for;assert a mod d eq 0;Append(~cl,a div d);
end for;
assert cl eq [ns[d]+#[x:x in ds|x eq d]:d in [1..5]];
Write("work/deg10/continuation_20261029/cont_local_6_12_q7_census.m",Sprint([*6,12,7,ds,ns,Coefficients(f),cl*]));
print "DECIC_LOCAL_CENSUS_COMPLETE",6,12,7,ns,Cputime(t0);
AF:=AlgorithmicFunctionField(FunctionField(C));c0:=[c:c in cs|Degree(c) eq 1][1];
G,_,phi:=ClassGroup(AF:ReductionDivisor:=FunctionFieldDivisor(Divisor(c0)),Proof:=true);
print "DECIC_LOCAL_CLASSGROUP_PROVEN",Invariants(G),Cputime(t0);
quit;
