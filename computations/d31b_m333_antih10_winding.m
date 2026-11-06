SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(3*10^9);tt:=Cputime();
Alarm(550);
for row in [[99,34,81]] do
 N,hh,genus:=Explode(row);D:=FullDirichletGroup(N);total:=0;packet:=[* *];index:=0;
 for chi in GaloisConjugacyRepresentatives(D) do
  if IsOdd(chi) or Evaluate(chi,hh) ne 1 or Evaluate(chi,10) ne -1 then continue;end if;index+:=1;
  S:=CuspidalSubspace(ModularSymbols(MinimalBaseRingCharacter(chi),2,1));d:=Dimension(S);total+:=d*EulerPhi(Order(chi));if d eq 0 then continue;end if;
  vambient:=Vector(Eltseq(WindingElement(S)))*ProjectionMatrix(S);v:=Vector(Coordinates(VectorSpace(S),vambient));
  ops:=[HeckeOperator(S,q):q in [2,3,5,7,11,13,17,19,23]];K:=BaseRing(ops[1]);rat:=K cmpeq Rationals();cs:=&cat[Eltseq(K!x):x in Eltseq(v) cat &cat[Eltseq(A):A in ops]];prime:=0;
  for mult in [1..5000] do
   p:=Order(chi)*mult+1;if p lt 5 or not IsPrime(p) or not &and[Denominator(c) mod p ne 0:c in cs] then continue;end if;F:=GF(p);
   if rat then r:=F!0;else rr:=Roots(PolynomialRing(F)!DefiningPolynomial(K));if #rr eq 0 then continue;end if;r:=rr[1][1];end if;
   prime:=p;break;
  end for;assert prime gt 0;
  red:=func<x|&+[F!c*r^(i-1):i->c in Eltseq(K!x)]>;
  vm:=Vector(F,[red(x):x in Eltseq(v)]);ams:=[Matrix(F,d,d,[red(x):x in Eltseq(A)]):A in ops];V:=VectorSpace(F,d);U:=sub<V|[vm]>;
  repeat old:=Dimension(U);U:=sub<V|Basis(U) cat [b*A:b in Basis(U),A in ams]>;until Dimension(U) eq old or Dimension(U) eq d;
  assert Dimension(U)eq d;print "D31B_M333_WINDING_COMPONENT",index,Order(chi),d,Dimension(U);Append(~packet,[*index,Order(chi),d,prime,[Integers()!x:x in Eltseq(vm)],[[Integers()!x:x in Eltseq(A)]:A in ams]*]);
  print "WINDING_COMPONENT_CHECKED",N,hh,index,Order(chi),d,prime,Cputime(tt);
 end for;
 assert total eq 40 and Dimension(JH(N,[hh])) eq genus and Dimension(JH(N,[hh,10]))eq 41;assert &and[Order(chi)in[2,10]:chi in GaloisConjugacyRepresentatives(D)|not IsOdd(chi)and Evaluate(chi,hh)eq 1 and Evaluate(chi,10)eq -1];
 Write(Sprintf("work/deg10/continuation_20261031b/M3N33/winding_antih10_N%o_H%o.m",N,hh),Sprint([*N,hh,genus,total,packet*]));
 print "WINDING_RANK_CERTIFICATE_COMPLETE",N,hh,genus,Cputime(tt);
end for;
print "D31B_M333_ANTIH10_WINDING_COMPLETE",Cputime(tt);quit;
