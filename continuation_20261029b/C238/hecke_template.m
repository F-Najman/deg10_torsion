SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(4*10^9);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");Attach("work/deg10/continuation_20261029/deps/fast_hecke.m");Attach("work/deg10/continuation_20261029b/deps/fast_hecke_cyclic_commonzero.m");
load "work/deg10/continuation_20261029b/deps/pair_af_place_linear.m";
Y:=MDX1(38,GF(3));C:=Curve(Y);assert Genus(C) eq 28;cs:=Cusps(Y);assert {*Degree(c):c in cs*} eq {*1^^18,9^^2*};FF:=FunctionField(C);AF,mp:=AlgorithmicFunctionField(FF);coords:=[mp(FF.i):i in [1,2]];
data:=eval Read("work/deg10/continuation_20261029b/C238/filter19_v2.m");source:=data[2];small:=data[3];remaining:=data[7];assert #remaining eq 35;
own:=AssociativeArray();sm:=AssociativeArray();fields:=AssociativeArray();R:=PolynomialRing(GF(3));
for d in [4..10] do own[d]:=[r:r in source|r[1] eq d];if d le 6 then sm[d]:=[r:r in small|r[1] eq d];end if;a:=data[4][d-3];assert a[1] eq d;fields[d]:=ext<GF(3)|R!a[2]>;end for;
cache:=AssociativeArray();images:=[* *];tests:=[* *];left:=[* *];
procedure Push(~DD,~cache,~images,row,mapno)
 d:=row[1];ky:=row[3];kk:=[mapno,d] cat ky;if IsDefined(cache,kk) then DD:=cache[kk];return;end if;
 K:=fields[d];z:=K.1;b:=Evaluate(R!ky[1..d],z);c:=Evaluate(R!ky[d+1..2*d],z);E:=EllipticCurve([K|0,c,0,b*(1-b-c),0]);P:=E![0,0];Q:=E![b,b];assert Order(Q) eq 38 and P ne 19*Q;
 if mapno le 2 then EG:=E;U:=mapno eq 1 select Q else P+Q;else T:=mapno eq 3 select P else P+19*Q;EG,ph:=IsogenyFromKernel(E,PolynomialRing(K).1-T[1]);U:=ph(Q);end if;assert Order(U) eq 38;
 xy:=TateNormalForm_xy(EG,U);pl:=C48PairAFPlaceLinear(coords,xy);e:=Degree(pl);assert d mod e eq 0;w:=d div e;cp:=CurvePlace(C,pl);DD:=w*Divisor(cp);cache[kk]:=DD;Append(~images,[*kk,e,w,[Eltseq(v):v in xy]*]);
end procedure;
for k in [1..#remaining] do
 if (k-1) mod 4 ne SHARD then continue;end if;desc:=remaining[k];sp:=[*own[desc[2]][desc[3]]*];if desc[1] eq 2 then Append(~sp,sm[desc[4]][desc[5]]);end if;cd:=10-&+[r[1]:r in sp];assert cd in [0..5];checks:=[* *];passed:=true;
 for mapno in [1..4] do
  N:=DivisorGroup(C)!0;for row in sp do DD:=DivisorGroup(C)!0;Push(~DD,~cache,~images,row,mapno);N+:=DD;end for;
  for q in [5,7] do
   T:=HeckeOperatorFast(Y,q,N);pp,ee:=Support(N);DN:=&+[ee[i]*Divisor(DiamondOperator(Y,q,pp[i])):i in [1..#pp]];A:=T-N-q*DN;assert Degree(A) eq 0 and Degree(T) eq (q+1)*(10-cd);
   yes:=IsPrincipal(A);Append(~checks,[*mapno,q,Degree(T),yes*]);print "C238_HECKE_TEST",k,desc,mapno,q,yes,Cputime(tt);if not yes then passed:=false;break;end if;
  end for;
  if not passed then break;end if;
 end for;
 Append(~tests,[*k,desc,cd,checks,passed*]);if passed then Append(~left,desc);end if;
 Write(Sprintf("work/deg10/continuation_20261029b/C238/hecke_shard%o.m",SHARD),Sprint([*images,tests,left*]):Overwrite:=true);
end for;
print "C238_HECKE_SHARD_COMPLETE",SHARD,#tests,#left,Cputime(tt);quit;
