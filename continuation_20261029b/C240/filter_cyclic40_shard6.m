/* Fresh degree9 seven-case test using only archived degree-independent
   exact Picard geometry, with new direct ideal identities for unknown primes. */
SetColumns(0);SetNthreads(1);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
load "work/deg10/continuation_20261029b/deps/pair_af_place_linear.m";
load "work/deg10/continuation_20261029b/C240/portable/portable_inputs.m";
load "work/deg10/continuation_20261029b/C240/portable/portable_forward.m";
k:=GF(3);PP:=PolynomialRing(k);B<t>:=FunctionField(k);RF<Y>:=PolynomialRing(B);
function Rat(v) return B!(PP!v[1])/B!(PP!v[2]);end function;
f:=RF![Rat(v):v in PortableFieldEquation];f/:=LeadingCoefficient(f);assert IsMonic(f);F<y>:=FunctionField(f);assert Genus(F) eq 25;
function GetFunction(data) return F![Rat(v):v in data];end function;
coords:=[GetFunction(v):v in PortableSutherland];tatecoords:=[GetFunction(v):v in PortableTate];
OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);fb:=[];
for j in [1..#PortableIdealBases] do
 OO:=PortableFiniteFlags[j] eq 1 select OF else OI;I:=ideal<OO|[GetFunction(v):v in PortableIdealBases[j]]>;pl:=Place(I);assert Degree(pl) eq PortableDegrees[j];Append(~fb,pl);
end for;assert #fb eq 616;print "C240_PORTABLE_PROVEN_GEOMETRY_READY",Cputime(t0);
Z:=FreeAbelianGroup(9);G,q:=quo<Z|[Z!v:v in PortableRelations]>;fbclasses:=[q(Z!v):v in PortableClasses];A:=Divisor(fb[609]);assert Degree(A) eq 1 and fbclasses[609] eq q(Z.9);
assert &*[v:v in Invariants(G)|v ne 0] eq 31579088486400;
QZ:=PolynomialRing(Integers());assert Evaluate(QZ!PortableFrobenius,1) eq 31579088486400;
phi,witness:=C240PortableDirectRR(F,fb,fbclasses,G,A,"work/deg10/continuation_20261029b/C240/portable/shard6_RR_identity.log");
hg:=[q(Z!v):v in PortableRationalCuspGenerators];H:=sub<G|hg>;assert Invariants(H) eq [2,2,2,8,120,895440];QH,qh:=quo<G|H>;
cc:=[q(Z!v):v in PortableCuspClasses];comps:=[*{QH!0}*] cat [*{QH|}:i in [1..6]*];for j in [1..#cc] do w:=PortableCuspDegrees[j];if w le 6 then for d in [w..6] do comps[d+1] join:={z+qh(cc[j]):z in comps[d-w+1]};end for;end if;end for;
groupdata:=[*Invariants(G),[Eltseq(r):r in Rows(RelationMatrix(G))],[Eltseq(v):v in hg],PortableCuspDegrees,[Eltseq(v):v in cc],[[Eltseq(z):z in co]:co in comps]*];
print "C240_PORTABLE_PROVEN_CLASS_MAP_READY",Invariants(G),Invariants(H),[#co:co in comps],Cputime(t0);
data:=eval Read("work/deg10/continuation_20261029b/C240/filter_quotient20_linear_arithmetic.txt");singles:=data[3];small:=data[4];fields:=data[2];
survivors:=(eval Read("work/deg10/continuation_20261029b/C240/filter_quotient20_isogeny_arithmetic.txt"))[4];
own:=AssociativeArray();sm:=AssociativeArray();for d in [4..10] do own[d]:=[r:r in singles|r[1] eq d];end for;for d in [4,5,6] do sm[d]:=[r:r in small|r[1] eq d];end for;
R:=PolynomialRing(k);fieldobjects:=AssociativeArray();for rec in fields do fieldobjects[rec[1]]:=ext<k|R!rec[2]>;end for;
cache:=AssociativeArray();images:=[* *];tests:=[* *];remaining:=[* *];
procedure Push(~answer,~cache,~images,row,mapno)
 d:=row[1];ky:=row[3];cachekey:=[mapno,d] cat ky;if IsDefined(cache,cachekey) then answer:=cache[cachekey];return;end if;
 K:=fieldobjects[d];z:=K.1;b:=Evaluate(R!ky[1..d],z);c:=Evaluate(R!ky[d+1..2*d],z);E:=EllipticCurve([K|0,c,0,b*(1-b-c),0]);P:=E![0,0];Q:=E![b,b];assert 40*Q eq E!0 and 20*Q ne E!0 and 8*Q ne E!0 and P ne 20*Q;
 if mapno le 2 then EG:=E;QG:=mapno eq 1 select Q else P+Q;
 else T:=mapno eq 3 select P else P+20*Q;EG,iso:=IsogenyFromKernel(E,PolynomialRing(K).1-T[1]);QG:=iso(Q);end if;
 assert 40*QG eq EG!0 and 20*QG ne EG!0 and 8*QG ne EG!0;
 xy:=TateNormalForm_xy(EG,QG);pl:=C48PairAFPlaceLinear(coords,xy);e:=Degree(pl);assert d mod e eq 0;w:=d div e;answer:=w*phi(Divisor(pl)-e*A);cache[cachekey]:=answer;
 Append(~images,[*cachekey,e,w,Eltseq(answer)*]);
end procedure;
for j in [1..#survivors] do
 if (j-1) mod 8 ne 6 then continue;end if;
 desc:=survivors[j];sp:=[*own[desc[2]][desc[3]]*];if desc[1] eq 2 then Append(~sp,sm[desc[4]][desc[5]]);end if;
 cd:=10-&+[r[1]:r in sp];assert cd in [0..6];passed:=true;checks:=[* *];
 for mapno in [1..4] do
  cl:=G!0;for row in sp do aa:=G!0;Push(~aa,~cache,~images,row,mapno);cl+:=aa;end for;
  yes:=-qh(cl) in comps[cd+1];Append(~checks,[*mapno,Eltseq(cl),yes*]);if not yes then passed:=false;break;end if;
 end for;
 Append(~tests,[*j,desc,cd,checks,passed*]);if passed then Append(~remaining,desc);end if;
 Write("work/deg10/continuation_20261029b/C240/portable/shard6_arithmetic.txt",Sprint([*groupdata,images,tests,remaining*]):Overwrite:=true);
 print "C240_PORTABLE_SEVEN_CASE",j,desc,passed,Cputime(t0);
end for;
fprintf witness,"]\n";Flush(witness);print "C240_D10_PORTABLE_CYCLIC40_SHARD6_DONE",#tests,#remaining,"complete=true",Cputime(t0);quit;
