SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(4*10^9);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
Attach("work/deg10/continuation_20261029/deps/fast_hecke.m");
Attach("work/deg10/shared/deps/class_rational_reduction.m");
load "work/deg10/continuation_20261029b/deps/pair_af_place_linear.m";
load "work/deg10/continuation_20261029b/deps/class_group_direct_RR_checked_v4_fastfb.m";
Y:=MDX1(25,GF(3));C:=Curve(Y);assert Genus(C) eq 12;os:=[Type(o) eq PlcCrvElt select Divisor(o) else o:o in CuspOrbitsQ(Y)];cs:=Cusps(Y);c0:=Support([o:o in os|Degree(o) eq 1][1])[1];
G,m1,m2:=ClassGroupRationalReduction(C,c0:Proof:=true);AF,mp:=AlgorithmicFunctionField(FunctionField(C));base:=FunctionFieldDivisor(Divisor(c0));coords:=[mp(FunctionField(C).i):i in [1,2]];
GG,_,phi,witness:=MDClassGroupDirectRRCheckedV4FastFB(AF,"work/deg10/continuation_20261029b/C250/target25_RR.log");assert GG cmpeq G;
hg:=[phi(FunctionFieldDivisor(o-Degree(o)*Divisor(c0))):o in os];H:=sub<G|hg>;print "C250_EXACT_CUSP_GROUP",Invariants(H),#H;QH,qh:=quo<G|H>;
cc:=[phi(FunctionFieldDivisor(Divisor(c)-Degree(c)*Divisor(c0))):c in cs];allowed:=[*{QH!0}*] cat [*{QH|}:i in [1..6]*];
for i in [1..#cs] do w:=Degree(cs[i]);if w gt 6 then continue;end if;for d in [w..6] do allowed[d+1] join:={z+qh(cc[i]):z in allowed[d-w+1]};end for;end for;
groupdata:=[*Invariants(G),[Eltseq(r):r in Rows(RelationMatrix(G))],[Eltseq(x):x in hg],[Degree(c):c in cs],[Eltseq(x):x in cc],[[Eltseq(x):x in a]:a in allowed]*];
print "C250_TARGET19_GROUP",Invariants(G),Invariants(H),[#a:a in allowed],Cputime(tt);
R:=PolynomialRing(GF(3));fields:=AssociativeArray();own:=AssociativeArray();sm:=AssociativeArray();source:=[* *];small:=[* *];coverage:=[* *];
function Key(v,D)
 ims:=[];for i in [1..D] do Append(~ims,&cat[[Integers()!a:a in Eltseq(z,GF(3))]:z in v]);v:=[z^3:z in v];end for;return Sort(ims)[1];
end function;
function BC(E,P,Q)
 ai:=aInvariants(E);xx:=Q[1]-P[1];yy:=Q[2];return [xx^3/yy^2,(ai[2]+3*P[1])*xx^2/yy^2];
end function;
for d in [6..10] do
 a:=eval Read(Sprintf("work/deg10/continuation_20261029b/C250/source%o.m",d));assert a[1] eq d and a[3];K<z>:=ext<GF(3)|R!a[2]>;fields[d]:=K;own[d]:=[* *];allkeys:={};cover:=[* *];
 for r in a[7] do
  if r[1] ne d then continue;end if;i:=#own[d]+1;ky:=r[2];row:=[*d,i,ky*];Append(~own[d],row);Append(~source,row);
 end for;
 Append(~coverage,[*d,a[2],cover*]);print "C250_FILTER_SOURCE",d,#own[d],Cputime(tt);
end for;
cache:=AssociativeArray();images:=[* *];tests:=[* *];remaining:=[* *];
procedure Push(~answer,~cache,~images,row,mapno)
 d:=row[1];ky:=row[3];kk:=[mapno,d] cat ky;if IsDefined(cache,kk) then answer:=cache[kk];return;end if;
 K:=fields[d];z:=K.1;b:=Evaluate(R!ky[1..d],z);c:=Evaluate(R!ky[d+1..2*d],z);E:=EllipticCurve([K|0,c,0,b*(1-b-c),0]);P:=E![0,0];Q:=E![b,b];assert Order(Q) eq 50 and P ne 25*Q;V:=2*Q;
 if mapno eq 1 then EG:=E;U:=V;else T:=[P,25*Q,P+25*Q][mapno-1];EG,ph:=IsogenyFromKernel(E,PolynomialRing(K).1-T[1]);U:=ph(V);end if;assert Order(U) eq 25;
 pair:=TateNormalForm_xy(EG,U);pl:=C48PairAFPlaceLinear(coords,pair);e:=Degree(pl);assert d mod e eq 0;w:=d div e;answer:=w*phi(Divisor(pl)-e*base);cache[kk]:=answer;Append(~images,[*kk,e,w,Eltseq(answer)*]);
end procedure;
supports:=[*[1,r[1],r[2]]:r in source*];
for k in [1..#supports] do
 desc:=supports[k];sp:=[*own[desc[2]][desc[3]]*];if desc[1] eq 2 then Append(~sp,sm[desc[4]][desc[5]]);end if;cd:=10-&+[r[1]:r in sp];assert cd in [0..6];passed:=true;checks:=[* *];
 for mapno in [1..4] do
  cl:=G!0;for r in sp do v:=G!0;Push(~v,~cache,~images,r,mapno);cl+:=v;end for;yes:=-qh(cl) in allowed[cd+1];Append(~checks,[*mapno,Eltseq(cl),yes*]);if not yes then passed:=false;break;end if;
 end for;
 Append(~tests,[*k,desc,cd,checks,passed*]);if passed then Append(~remaining,desc);end if;
 if k mod 50 eq 0 then print "C250_TARGET19_PROGRESS",k,#remaining,Cputime(tt);end if;
end for;
fprintf witness,"]\n";Flush(witness);
Write("work/deg10/continuation_20261029b/C250/filter25.m",Sprint([*groupdata,source,small,coverage,images,tests,remaining*]));
print "C250_TARGET19_COMPLETE",#supports,#remaining,Cputime(tt);quit;
