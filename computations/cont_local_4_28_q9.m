SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(5*10^9);t0:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
Attach("work/deg10/continuation_20261029/deps/fast_hecke.m");
Attach("work/deg10/continuation_20261029/deps/fast_hecke_11.m");
X:=MDX11(4,28,GF(9));
C:=Curve(X);assert Genus(C) eq 61;
cs:=Cusps(X);ds:=[Degree(c):c in cs];
print "DECIC_LOCAL_CUSPS",4,28,9,ds,Cputime(t0);
ns:=[];
for d in [1..5] do xs:=NoncuspidalPlaces(X,d);Append(~ns,#xs);print "DECIC_LOCAL_COUNT",4,28,9,d,#xs,Cputime(t0);end for;
J:=JH(112,[29]);assert Dimension(J) eq 61;
p:=Characteristic(GF(9));e:=Degree(GF(9));f:=FrobeniusPolynomial(J,p);assert Degree(f) eq 122;tr:=[];
for i in [1..5*e] do
 a:=-i*Coefficient(f,122-i);
 for j in [1..i-1] do a-:=Coefficient(f,122-j)*tr[i-j];end for;
 Append(~tr,Integers()!a);
end for;
cl:=[];
for d in [1..5] do
 a:=9^d+1-tr[e*d];for j in Divisors(d) do if j lt d then a-:=j*cl[j];end if;end for;assert a mod d eq 0;Append(~cl,a div d);
end for;
assert cl eq [ns[d]+#[x:x in ds|x eq d]:d in [1..5]];
Write("work/deg10/continuation_20261029/cont_local_4_28_q9_census.m",Sprint([*4,28,9,ds,ns,Coefficients(f),cl*]));
print "DECIC_LOCAL_CENSUS_COMPLETE",4,28,9,ns,Cputime(t0);
quit;
