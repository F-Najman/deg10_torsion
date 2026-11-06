SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(4*10^9);tt:=Cputime();
load "work/deg10/continuation_20261029b/deps/pair_af_place_linear.m";
load "work/deg10/continuation_20261029b/C240/portable/portable_forward.m";
D:=eval Read("work/deg10/continuation_20261029b/C226/proven_picard.m");k:=GF(3);PP:=PolynomialRing(k);B<t>:=FunctionField(k);RF<Y>:=PolynomialRing(B);
function Rat(v) return B!(PP!v[1])/B!(PP!v[2]);end function;
f:=RF![Rat(v):v in D[1]];f/:=LeadingCoefficient(f);F<y>:=FunctionField(f);assert Genus(F) eq 25;
function GetFunction(data) return F![Rat(v):v in data];end function;
OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);fb:=[];
for j in [1..#D[4]] do OO:=D[2][j] eq 1 select OF else OI;I:=ideal<OO|[GetFunction(v):v in D[4][j]]>;pl:=Place(I);assert Degree(pl) eq D[3][j];Append(~fb,pl);end for;
Z:=FreeAbelianGroup(#D[5][1]);G,q:=quo<Z|[Z!v:v in D[6]]>;fc:=[q(Z!v):v in D[5]];A:=Divisor(fb[D[8]]);coords:=[GetFunction(v):v in D[7]];
hg:=[fc[j]-fc[D[8]]:j in D[9]];cc:=[fc[j]-3*fc[D[8]]:j in D[10]];H:=sub<G|hg>;assert #H eq D[12] and Invariants(H) eq D[13];QH,qh:=quo<G|H>;
phi,witness:=C240PortableDirectRR(F,fb,fc,G,A,Sprintf("work/deg10/continuation_20261029b/C226/picard%o_RR.log",SHARD));
group:=[*Invariants(G),[Eltseq(v):v in Rows(RelationMatrix(G))],[Eltseq(v):v in hg],[Eltseq(v):v in cc]*];
data:=eval Read("work/deg10/continuation_20261029b/C226/filter26.m");source:=data[2];small:=data[3];remaining:=data[7];assert #remaining eq 32;own:=AssociativeArray();sm:=AssociativeArray();fields:=AssociativeArray();
for d in [5..10] do own[d]:=[r:r in source|r[1] eq d];if d le 6 then sm[d]:=[r:r in small|r[1] eq d];end if;a:=data[4][d-4];assert a[1] eq d;fields[d]:=ext<k|PP!a[2]>;end for;
cache:=AssociativeArray();images:=[* *];tests:=[* *];left:=[* *];
procedure Push(~cl,~cache,~images,row)
 d:=row[1];ky:=row[3];kk:=[d] cat ky;if IsDefined(cache,kk) then cl:=cache[kk];return;end if;K:=fields[d];z:=K.1;bc:=[Evaluate(PP!ky[1..d],z),Evaluate(PP!ky[d+1..2*d],z)];pl:=C48PairAFPlaceLinear(coords,bc);assert Degree(pl) eq d;cl:=phi(Divisor(pl)-d*A);cache[kk]:=cl;Append(~images,[*kk,Eltseq(cl)*]);
end procedure;
for i in [1..#remaining] do
 if (i-1) mod 8 ne SHARD then continue;end if;desc:=remaining[i];rows:=[*own[desc[2]][desc[3]]*];if desc[1] eq 2 then Append(~rows,sm[desc[4]][desc[5]]);end if;cd:=10-&+[r[1]:r in rows];cl:=G!0;
 for row in rows do v:=G!0;Push(~v,~cache,~images,row);cl+:=v;end for;
 allowed:={QH!0};if cd ge 3 then allowed join:={-qh(c):c in cc};end if;yes:=qh(cl) in allowed;Append(~tests,[*i,desc,cd,Eltseq(cl),yes*]);if yes then Append(~left,desc);end if;
 Write(Sprintf("work/deg10/continuation_20261029b/C226/picard%o.m",SHARD),Sprint([*group,images,tests,left*]):Overwrite:=true);print "C226_PICARD_TEST",i,desc,yes,Cputime(tt);
end for;
fprintf witness,"]\n";Flush(witness);print "C226_PICARD_SHARD_COMPLETE",SHARD,#tests,#left,Cputime(tt);quit;
