SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(4*10^9);tt:=Cputime();AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
load "work/deg10/next_cases_20261030/C242/geometry/deps/portable_forward.m";
load "work/deg10/next_cases_20261030/C242/geometry/quotient5_helpers.m";
a:=eval Read("work/deg10/next_cases_20261030/C242/geometry/x0_42_proven_picard_v2.m");k:=GF(5);PP:=PolynomialRing(k);B<t>:=FunctionField(k);RF<Y>:=PolynomialRing(B);
function Rat(v) return B!(PP!v[1])/B!(PP!v[2]);end function;
f:=RF![Rat(v):v in a[1]];f/:=LeadingCoefficient(f);F<y>:=FunctionField(f);assert Genus(F) eq 5;
function GFUN(v) return F![Rat(w):w in v];end function;
OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
function DP(s) OO:=s[1] eq 1 select OF else OI;I:=ideal<OO|[GFUN(v):v in s[3]]>;pl:=Place(I);assert Degree(pl) eq s[2];return pl;end function;
fb:=[];for j in [1..#a[4]] do OO:=a[2][j] eq 1 select OF else OI;I:=ideal<OO|[GFUN(v):v in a[4][j]]>;pl:=Place(I);assert Degree(pl) eq a[3][j];Append(~fb,pl);end for;
Zn:=FreeAbelianGroup(#a[6][1]);G,q:=quo<Zn|[Zn!v:v in a[6]]>;fc:=[q(Zn!v):v in a[5]];A0:=Divisor(fb[a[8]]);A:=Divisor(DP(a[12]));coords:=[GFUN(v):v in a[7]];H:=sub<G|[q(Zn!v):v in a[13]]>;assert Invariants(G) eq [2,4,48,48,0] and Invariants(H) eq [2,48];QH,qh:=quo<G|H>;assert #TorsionSubgroup(QH) eq 192;
cc:=[q(Zn!v):v in a[10]];comps:=[*{QH!v:v in co}:co in a[14]*];recomp:=[*{QH!0}*] cat [*{QH|}:d in [1..6]*];
for i in [1..#cc] do w:=a[9][i];if w le 6 then for d in [w..6] do recomp[d+1] join:={z+qh(cc[i]):z in recomp[d-w+1]};end for;end if;end for;assert comps eq recomp and [#s:s in comps] eq [1,6,12,16,20,24,24];
phi,witness:=C240PortableDirectRR(F,fb,fc,G,A0,OUTPUTPREFIX cat "_RR.log");
data:=eval Read(INPUTPATH);fields:=AssociativeArray();for row in data[1] do fields[row[1]]:=ext<k|PP!row[2]>;end for;cases:=data[2];
cache:=AssociativeArray();images:=[* *];tests:=[* *];remaining:=[* *];
function NormalizePoint(E,S)
 ai:=aInvariants(E);assert ai[1] eq 0 and ai[3] eq 0;x:=S[1];y:=S[2];m:=(3*x^2+2*ai[2]*x+ai[4])/(2*y);a2p:=3*x+ai[2]-m^2;a3p:=2*y;kap:=a3p/a2p;bb:=-a2p^3/a3p^2;cc:=1-2*m*a2p/a3p;assert [bb,cc] eq TateNormalForm_bc(E,S);return bb,cc,kap;
end function;
function ProjectX0(E,Q)
 b6,c6,k6:=NormalizePoint(E,7*Q);b7,c7,k7:=NormalizePoint(E,6*Q);assert b6 eq c6+c6^2;v:=b7/c7;assert b7 eq v^2*(v-1) and c7 eq v*(v-1);t:=(v^3-3*v+1)/(v*(v-1))-8;H:=v^2-v+1;rat:=(k7*H/k6)^2;D:=t^2+13*t+49;N:=t^2+5*t+1;M:=t^4+14*t^3+63*t^2+70*t-7;EE6:=EllipticCurve([BaseRing(E)|1-c6,-b6,-b6,0,0]);assert jInvariant(E) eq N^3*D/t and rat*cInvariants(EE6)[1]*M+cInvariants(EE6)[2]*N*D eq 0;return[c6,t,rat];
end function;
procedure ProjectClass(~answer,~cache,~images,row,mapno)
 d:=row[1];ky:=row[2];ck:=[mapno,d] cat ky;if IsDefined(cache,ck) then answer:=cache[ck];return;end if;K:=fields[d];z:=K.1;b:=Evaluate(PP!ky[1..d],z);c:=Evaluate(PP!ky[d+1..2*d],z);E:=EllipticCurve([K|0,c,0,b*(1-b-c),0]);P:=E![0,0];Q:=E![b,b];assert Order(Q) eq 42 and P ne 21*Q;
 if mapno le 3 then EG:=E;QG:=mapno eq 1 select Q else(mapno eq 2 select P+Q else P+22*Q);
 else T:=mapno eq 4 select P else(mapno eq 5 select P+21*Q else 21*Q);EG,iso:=IsogenyFromKernel(E,PolynomialRing(K).1-T[1]);QG:=mapno eq 6 select iso(P)+iso(Q) else iso(Q);end if;assert Order(QG) eq 42;
 vals:=ProjectX0(EG,QG);pl:=InverseThree(coords,vals);e:=Degree(pl);assert d mod e eq 0;mult:=d div e;answer:=mult*phi(Divisor(pl)-e*A);cache[ck]:=answer;Append(~images,[*mapno,d,ky,e,mult,Eltseq(answer),[[Integers()!x:x in Eltseq(v,k)]:v in vals]*]);
end procedure;
limit:=Minimum(#cases,MAXCASES);
for index in [1..limit] do
 row:=cases[index];id:=row[1];points:=row[2];ndeg:=&+[p[1]:p in points];cd:=10-ndeg;assert cd in [0..6];passed:=true;checks:=[* *];
 for mapno in [1..6] do
  cl:=G!0;for p in points do a0:=G!0;ProjectClass(~a0,~cache,~images,p,mapno);cl+:=a0;end for;
  yes:=-qh(cl) in comps[cd+1];Append(~checks,[*mapno,Eltseq(cl),yes*]);if not yes then passed:=false;break;end if;
 end for;
 Append(~tests,[*id,cd,checks,passed*]);if passed then Append(~remaining,row);end if;
 if index mod 10 eq 0 or index eq limit then
  Write(OUTPUTPREFIX cat "_arithmetic.m",Sprint([*data[1],images,tests,remaining*]):Overwrite:=true);
  printf "C242_D10_X0_42_STAGE3_PROGRESS tested=%o total=%o remaining=%o imagecount=%o cpu=%o\n",index,limit,#remaining,#images,Cputime(tt);
 end if;
end for;
fprintf witness,"]\n";Flush(witness);printf "C242_D10_X0_42_STAGE3_COMPLETE tested=%o total=%o remaining=%o images=%o full_input=%o cpu=%o\n",limit,#cases,#remaining,#images,limit eq #cases,Cputime(tt);quit;
