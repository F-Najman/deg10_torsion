/* Exact helpers from the accepted genus-nine quotient certificate.
   Native coordinate order is [b12,c12,w]. Requires no source changes. */
function KeyThree(vs,K,p,d)
 keys:=[];for j in [1..d] do Append(~keys,&cat[[Integers()!a:a in Eltseq(v,GF(p))]:v in vs]);vs:=[v^p:v in vs];end for;return Sort(keys)[1];
end function;
function InverseThree(coords,values)
 F:=Parent(coords[1]);K:=Parent(values[1]);k:=PrimeField(K);p:=#k;e:=LCM([Degree(MinimalPolynomial(v,k)):v in values]);found:=false;
 for bits in [[1,0,0],[0,1,0],[0,0,1]] cat [[1,a,b]:a,b in [0..p-1]] cat [[0,1,a]:a in [0..p-1]] do
  gen:=&+[K!bits[i]*values[i]:i in [1..3]];pol:=MinimalPolynomial(gen,k);if Degree(pol) eq e then chosen:=bits;found:=true;break;end if;
 end for;assert found;M:=Matrix(k,[Eltseq(gen^i,k):i in [0..e-1]]);assert Rank(M) eq e;
 cc:=[Eltseq(Solution(M,Vector(k,Eltseq(v,k)))):v in values];assert &and[&+[K!cc[j][i+1]*gen^i:i in [0..e-1]] eq values[j]:j in [1..3]];
 gg:=&+[F!chosen[i]*coords[i]:i in [1..3]];fs:=[Evaluate(pol,gg)] cat [coords[j]-&+[F!cc[j][i+1]*gg^i:i in [0..e-1]]:j in [1..3]];
 O:=MaximalOrderFinite(F);I:=ideal<O|fs>;fa:=Factorization(I);places:=[Place(t[1]):t in fa|t[2] gt 0];
 if #places eq 0 then Oi:=MaximalOrderInfinite(F);Ii:=ideal<Oi|fs>;fi:=Factorization(Ii);places:=[Place(t[1]):t in fi|t[2] gt 0];end if;
 assert #places eq 1 and Degree(places[1]) eq e;assert &and[f eq 0 or Valuation(f,places[1]) gt 0:f in fs];return places[1];
end function;
function M2ThirtySixPair(E,P,Q)
 aa:=aInvariants(E);xx:=Q[1]-P[1];yy:=Q[2];return [xx^3/yy^2,(aa[2]+3*P[1])*xx^2/yy^2];
end function;
function Projection(E,P,Q)
 Q12:=3*Q;P3:=4*Q12;lambda:=Q12[2]/(Q12[1]-P[1]);bc:=M2ThirtySixPair(E,P,Q12);
 EG,iso:=IsogenyFromKernel(E,PolynomialRing(BaseRing(E)).1-P3[1]);S:=iso(4*Q);w:=(S[1]-P[1])/lambda^2;return bc cat [w];
end function;
