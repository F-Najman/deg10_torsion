SetColumns(0);SetQuitOnError(true);SetNthreads(1);Alarm(180);tt:=Cputime();AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");load "work/deg10/continuation_20261031b/Y228_shared/fast_primitive_inverse.m";load "work/deg10/continuation_20261031b/Y228_shared/union_cases.m";k:=GF(3);R<X>:=PolynomialRing(k);f:=R![2,2,-1,-2,-1,1];print "KNOWN_DECIC_QUINTIC_AT3",Factorization(f);good:=[* *];bad:=[];
for row in Factorization(f) do
 g,m:=Explode(row);d:=Degree(g);K<a>:=ext<k|g>;b:=Evaluate(R![1032,-239,-222,-759,419],a);c:=Evaluate(R![42,-11,-8,-31,17],a);ai:=[K|1-c,-b,-b,0,0];a1,a2,a3,a4,a6:=Explode(ai);b2:=a1^2+4*a2;b4:=a1*a3+2*a4;b6:=a3^2+4*a6;b8:=a1^2*a6+4*a2*a6-a1*a3*a4+a2*a3^2-a4^2;disc:=-b2^2*b8-8*b4^3-27*b6^2+9*b2*b4*b6;
 if disc eq 0 then Append(~bad,[d,m]);continue;end if;
 E:=EllipticCurve(ai);P:=E![0,0];assert Order(P) eq 28;T:=14*P;RR<Z>:=PolynomialRing(K);h:=DivisionPolynomial(E,2) div (Z-T[1]);assert Degree(h) eq 2 and IsSquarefree(h);
 for hh in Factorization(h) do pol,em:=Explode(hh);assert em eq 1;e:=Degree(pol);if e eq 1 then L:=K;r:=Roots(pol)[1][1];else L<r>:=ext<K|pol>;end if;EL:=ChangeRing(E,L);PL:=EL!P;S:=EL![r,-(L!a1*r+L!a3)/2];assert Order(S) eq 2 and S ne 14*PL;Append(~good,[*d,e,m,Y228PrimitiveJointKey(DerickxNormalForm_bc(EL,S,PL)),[Integers()!z:z in Coefficients(DefiningPolynomial(L,k))],[[Integers()!z:z in Eltseq(x,k)]:x in aInvariants(EL)],[[Integers()!z:z in Eltseq(x,k)]:x in[S[1],S[2],PL[1],PL[2]]]*]);
 end for;
end for;
print "KNOWN_DECIC_GOOD_BAD",[[z[1],z[2],z[3]]:z in good],bad;matches:=[* *];
for ci in[1..#UnionCases] do
 id,cusp,parts:=Explode(UnionCases[ci]);if #good ne #parts or cusp ne &+[2*z[1]*z[2]:z in bad] then continue;end if;
 targets:=[* *];for part in parts do dd,m,ky:=Explode(part);pp:=[z[2]:z in UnionFields|z[1] eq dd][1];KK<aa>:=ext<k|R!pp>;Append(~targets,[*dd,m,Y228PrimitiveJointKey([Evaluate(R!ky[1..dd],aa),Evaluate(R!ky[dd+1..2*dd],aa)])*]);end for;
 for bb in[0,1] do for cc in[0,1] do for u in[1,3,5,9,11,13] do
  keys:=[* *];for gr in good do L<aa>:=ext<k|R!gr[5]>;ev:=func<v|Evaluate(R!v,aa)>;EE:=EllipticCurve([ev(v):v in gr[6]]);ss:=EE![ev(gr[7][1]),ev(gr[7][2])];qq:=EE![ev(gr[7][3]),ev(gr[7][4])];Append(~keys,[*gr[1]*gr[2],gr[3],Y228PrimitiveJointKey(DerickxNormalForm_bc(EE,ss+14*bb*qq,cc*ss+u*qq))*]);end for;
  if Sort([Sprint(z):z in keys]) eq Sort([Sprint(z):z in targets]) then Append(~matches,[ci,bb,cc,u]);end if;
 end for;end for;end for;
end for;
Write("work/deg10/continuation_20261031b/Y228_shared/known_decic_reduction_at3.m",Sprint([*good,bad,matches*]):Overwrite:=true);print "D31B_Y228_KNOWN_DECIC_REDUCTION_COMPLETE",matches,Cputime(tt);quit;
