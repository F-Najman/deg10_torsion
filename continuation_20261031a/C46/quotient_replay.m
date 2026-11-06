SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(10*10^9);tt:=Cputime();
AttachSpec("work/degree9/shared/mdmagma/v2/mdmagma.spec");Attach("work/degree8/cyclic_sieves/class_rational_reduction.m");
load "work/deg10/continuation_20261031a/C46/all_source_cases.m";
k:=GF(3);PF:=PolynomialRing(k);X:=MDX1(46,k);Y:=MDX1(23,k);cf:=DefiningPolynomial(Curve(X));
cy:=Cusps(Y);c0:=[P:P in cy|Degree(P) eq 1][1];G,ig,cg:=ClassGroupRationalReduction(Curve(Y),c0:Proof:=true);
HH:=sub<G|[cg(Divisor(P)-Degree(P)*Divisor(c0)):P in cy]>;Q,qm:=quo<G|HH>;
assert Invariants(G) eq [9406793,0] and #HH eq 408991 and Invariants(Q) eq [23,0];
cache:=AssociativeArray();points:=[* *];caseresults:=[* *];rf:=recformat<P:PtEll>;
for a in AllSourceCases do
 typ,idx,oldreject,data:=Explode(a);total:=Q!0;refs:=[];
 for sr in data do
  d,mult,cp,co:=Explode(sr);ky:=Sprint([*d,cp,co*]);
  if not IsDefined(cache,ky) then
   L<ar>:=ext<k|PF!cp>;assert Degree(L) eq d;vv:=[Evaluate(PF!c,ar):c in co];assert Evaluate(cf,vv) eq 0 and vv[3] ne 0;
   x:=vv[1]/vv[3];y:=vv[2]/vv[3];assert Min([h:h in Divisors(d)|x^(3^h) eq x and y^(3^h) eq y]) eq d;
   rr:=(x^2*y-x*y+y-1)/(x^2*y-x);ss:=(x*y-y+1)/(x*y);b:=rr*ss*(rr-1);c:=ss*(rr-1);
   E:=EllipticCurve([L|1-c,-b,-b,0,0]);P:=E![0,0];assert Discriminant(E) ne 0 and 46*P eq E!0 and 23*P ne E!0 and 2*P ne E!0;
   P23:=2*P;assert 23*P23 eq E!0;ls:=rec<rf|P:=P23>;R:=ModuliPoint(Y,E,ls);dr:=Degree(R);assert d mod dr eq 0;
   cc:=cg(Divisor(R)-dr*Divisor(c0));v:=(d div dr)*qm(cc);assert Eltseq(v)[2] eq 0;
   xy:=TateNormalForm_xy(E,P23);row:=[*d,cp,co,[Eltseq(e,k):e in aInvariants(E)],dr,d div dr,[Eltseq(e,k):e in xy],Eltseq(cc),Eltseq(v)*];
   Append(~points,row);cache[ky]:=<#points,v>;
  end if;
  z:=cache[ky];total+:=mult*z[2];Append(~refs,<z[1],mult>);
 end for;
 reject:=total ne Q!0;assert reject eq oldreject;
 Append(~caseresults,[*typ,idx,refs,Eltseq(total),reject*]);
 if #caseresults mod 100 eq 0 then printf "D31A_C46_QUOTIENT_REPLAY_PROGRESS cases=%o points=%o time=%o\n",#caseresults,#points,Cputime(tt);end if;
end for;
assert #caseresults eq 1023 and #[a:a in caseresults|not a[5]] eq 37;
Write("work/deg10/continuation_20261031a/C46/quotient_replay_points.m",Sprintf("SourceTargetPoints:=%o; CaseClasses:=%o;",points,caseresults));
printf "D31A_C46_QUOTIENT_REPLAY_COMPLETE cases=%o points=%o quotient_rejects=986 exceptions=37 time=%o\n",#caseresults,#points,Cputime(tt);quit;
