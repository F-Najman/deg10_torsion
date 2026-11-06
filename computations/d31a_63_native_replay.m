SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(12*10^9);Alarm(800);tt:=Cputime();root:="work/deg10/continuation_20261031a/C63/";
k:=GF(2);PP:=PolynomialRing(k);K<t>:=FunctionField(k);R<Z>:=PolynomialRing(K);
function DR(v)return K!(PP!v[1])/K!(PP!v[2]);end function;
models:=[eval Read(root cat "good7_model.m"),eval Read(root cat "good10_model.m")];assert models[1][1] eq models[2][1] and models[1][2] eq models[2][2];
b:=t^2*(t-1)*(t^2-t+1);c:=t^2*(t-1);E9:=EllipticCurve([1-c,-b,-b,K!0,K!0]);f:=R![DR(q):q in models[1][1]];assert f eq R!DivisionPolynomial(E9,7) and Degree(f) eq 24 and IsIrreducible(f);F<u>:=FunctionField(f);assert Genus(F) eq 49;
function DF(v)return F![DR(q):q in v];end function;
jaf:=DF(models[1][2]);assert jaf eq F!jInvariant(E9);OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);Dzero:=DivisorGroup(F)!0;
function Key2(vals)
 L:=Parent(vals[1]);d:=LCM([Degree(MinimalPolynomial(v,k)):v in vals]);found:=false;
 for bits in [[1,0],[0,1],[1,1]] do gv:=&+[L!bits[i]*vals[i]:i in[1..2]];pol:=MinimalPolynomial(gv,k);if Degree(pol) eq d then chosen:=bits;found:=true;break;end if;end for;assert found;
 M:=Matrix(k,[Eltseq(gv^i,k):i in[0..d-1]]);assert Rank(M) eq d;co:=[Eltseq(Solution(M,Vector(k,Eltseq(v,k)))):v in vals];assert &and[&+[L!co[j][i+1]*gv^i:i in[0..d-1]] eq vals[j]:j in[1..2]];
 return [*d,chosen,[Integers()!q:q in Coefficients(pol)],[[Integers()!q:q in row]:row in co]*];
end function;
function Prime(pd)
 O:=pd[2] eq 1 select OF else OI;gs:=[DF(q):q in pd[3]];assert #gs eq 2;I:=ideal<O|gs>;assert IsIntegral(I) and IsPrime(I);Q:=Place(I);assert Ideal(Q) eq I and Degree(Q) eq pd[1];a,b:=TwoGenerators(Q);assert ideal<O|a,b> eq I;
 return Q;
end function;
function DecodeSD(rows)
 D:=Dzero;for row in rows do D+:=row[2]*Divisor(Prime(row[1]));end for;return D;
end function;
out:=[* *];total:=0;
for mi in[1,2] do
 model:=models[mi];e:=model[3];raw:=eval Read(root cat Sprintf("good%o_RR.m",e));assert #raw eq 3;L<z>:=ext<k|PP!model[5]>;cases:=[* *];
 for row in raw do
  ci:=row[1];P:=Prime(row[2]);Da:=Prime(row[3]);assert Degree(P) eq e and Degree(Da) eq e and Valuation(jaf,P) eq 0 and Valuation(jaf,Da) eq 0;
  keyP:=Key2([Evaluate(F!t,P),Evaluate(u,P)]);keyD:=Key2([Evaluate(F!t,Da),Evaluate(u,Da)]);
  T:=Dzero;hec:=[* *];assert &+[r[1]:r in row[4]] eq 6;
  for h in row[4] do Q:=Prime(h[3]);assert e*h[1] eq h[2]*Degree(Q);T+:=h[2]*Divisor(Q);Append(~hec,[*h[1],h[2],Key2([Evaluate(F!t,Q),Evaluate(u,Q)])*]);end for;
  assert Degree(T) eq 6*e;A:=T-Divisor(P)-5*Divisor(Da);h0s:=[];
  if e eq 10 then
   B:=DecodeSD(row[5]);assert B eq A;h0:=#Basis(B:Reduction:=false,Simplification:="None");assert h0 eq 0 and row[6] eq h0;Append(~h0s,h0);total+:=1;
  else
   assert e eq 7 and #row[5] eq 6;pairs:=[];
   for rr in row[5] do C1:=Prime(rr[3]);C2:=Prime(rr[4]);assert Degree(C1) eq 3 and Degree(C2) eq 3 and C1 ne C2 and Valuation(jaf,C1) eq -1 and Valuation(jaf,C2) eq -1;B:=A+4*Divisor(C1)-4*Divisor(C2);assert B eq DecodeSD(rr[5]);h0:=#Basis(B:Reduction:=false,Simplification:="None");assert h0 eq 0 and h0 eq rr[6];Append(~h0s,h0);Append(~pairs,<rr[1],rr[2]>);total+:=1;end for;
   assert Seqset(pairs) eq {<i,j>:i,j in[1..3]|i ne j};
  end if;
  Append(~cases,[*ci,keyP,keyD,hec,h0s*]);printf "D31A_C63_NATIVE_CASE e=%o ci=%o tests=%o time=%o\n",e,ci,#h0s,Cputime(tt);
 end for;
 Append(~out,[*e,cases*]);
end for;
assert total eq 21;Write(root cat "native_replay.m",Sprint(out):Overwrite:=true);printf "D31A_C63_NATIVE_REPLAY_COMPLETE tests=%o h0_all_zero=true two_generator_roundtrips=true time=%o\n",total,Cputime(tt);quit;
