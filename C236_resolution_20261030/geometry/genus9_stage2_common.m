SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(4*10^9);tt:=Cputime();AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
load "work/deg10/C236_resolution_20261030/geometry/deps/portable_forward.m";
load "work/deg10/C236_resolution_20261030/geometry/quotient9_projection_helpers.m";
a:=eval Read("work/deg10/C236_resolution_20261030/geometry/quotient9_proven_picard.m");k:=GF(5);PP:=PolynomialRing(k);B<t>:=FunctionField(k);RF<Y>:=PolynomialRing(B);
function Rat(v) return B!(PP!v[1])/B!(PP!v[2]);end function;
f:=RF![Rat(v):v in a[1]];f/:=LeadingCoefficient(f);F<y>:=FunctionField(f);assert Genus(F) eq 9;
function GFUN(v) return F![Rat(w):w in v];end function;
OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
function DP(s) OO:=s[1] eq 1 select OF else OI;I:=ideal<OO|[GFUN(v):v in s[3]]>;pl:=Place(I);assert Degree(pl) eq s[2];return pl;end function;
fb:=[];for j in [1..#a[4]] do OO:=a[2][j] eq 1 select OF else OI;I:=ideal<OO|[GFUN(v):v in a[4][j]]>;pl:=Place(I);assert Degree(pl) eq a[3][j];Append(~fb,pl);end for;
Zn:=FreeAbelianGroup(#a[6][1]);G,q:=quo<Zn|[Zn!v:v in a[6]]>;fc:=[q(Zn!v):v in a[5]];A0:=Divisor(fb[a[8]]);A:=Divisor(DP(a[12]));coords:=[GFUN(v):v in a[7]];H:=sub<G|[q(Zn!v):v in a[13]]>;assert Invariants(G) eq [3,3,3,24,24,192,0] and Invariants(H) eq [3,3,3,24,24,24];QH,qh:=quo<G|H>;assert #TorsionSubgroup(QH) eq 8;
phi,witness:=C240PortableDirectRR(F,fb,fc,G,A0,OUTPUTPREFIX cat "_RR.log");
data:=eval Read(INPUTPATH);fields:=AssociativeArray();for row in data[1] do fields[row[1]]:=ext<k|PP!row[2]>;end for;cases:=data[2];
cache:=AssociativeArray();images:=[* *];tests:=[* *];remaining:=[* *];
procedure ProjectClass(~answer,~cache,~images,row)
 d:=row[1];ky:=row[2];ck:=[d] cat ky;if IsDefined(cache,ck) then answer:=cache[ck];return;end if;
 K:=fields[d];z:=K.1;b:=Evaluate(PP!ky[1..d],z);c:=Evaluate(PP!ky[d+1..2*d],z);E:=EllipticCurve([K|0,c,0,b*(1-b-c),0]);P:=E![0,0];Q:=E![b,b];assert Order(Q) eq 36 and P ne 18*Q;
 v:=Projection(E,P,Q);pl:=InverseThree(coords,v);e:=Degree(pl);assert d mod e eq 0;mult:=d div e;assert mult in [1,3];answer:=mult*phi(Divisor(pl)-e*A);cache[ck]:=answer;
 Append(~images,[*d,ky,e,mult,Eltseq(answer),KeyThree(v,K,5,d)*]);
end procedure;
limit:=Minimum(#cases,MAXCASES);
for index in [1..limit] do
 row:=cases[index];id:=row[1];points:=row[2];ndeg:=&+[p[1]:p in points];assert ndeg in [3..10];cl:=G!0;
 for p in points do a0:=G!0;ProjectClass(~a0,~cache,~images,p);cl+:=a0;end for;
 yes:=cl in H;Append(~tests,[*id,Eltseq(cl),yes*]);if yes then Append(~remaining,row);end if;
 if index mod 20 eq 0 or index eq limit then
  Write(OUTPUTPREFIX cat "_arithmetic.m",Sprint([*data[1],images,tests,remaining*]):Overwrite:=true);
  printf "C236_D10_GENUS9_STAGE2_PROGRESS tested=%o total=%o remaining=%o imagecount=%o cpu=%o\n",index,limit,#remaining,#images,Cputime(tt);
 end if;
end for;
fprintf witness,"]\n";Flush(witness);printf "C236_D10_GENUS9_STAGE2_COMPLETE tested=%o total=%o remaining=%o images=%o full_input=%o cpu=%o\n",limit,#cases,#remaining,#images,limit eq #cases,Cputime(tt);quit;
