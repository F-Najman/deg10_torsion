/* Recover the small plane coordinates from the marked Derickx normal form,
   then construct the closed prime by a triangular finite-field ideal.
   Every returned place is bound to the complete physical (B,C) moduli key. */
function Y228PrimitiveJointKey(values)
 K:=Parent(values[1]);k:=PrimeField(K);p:=#k;e:=LCM([Degree(MinimalPolynomial(v,k)):v in values]);found:=false;
 for bits in [[1,0],[0,1]] cat [[1,a]:a in[1..p-1]] do
  gen:=&+[K!bits[i]*values[i]:i in[1,2]];pol:=MinimalPolynomial(gen,k);if Degree(pol) eq e then chosen:=bits;found:=true;break;end if;
 end for;assert found;M:=Matrix(k,[Eltseq(gen^i,k):i in[0..e-1]]);assert Rank(M) eq e;cs:=[Eltseq(Solution(M,Vector(k,Eltseq(v,k)))):v in values];
 return[*e,chosen,[Integers()!a:a in Coefficients(pol)],[[Integers()!a:a in row]:row in cs]*];
end function;
function Y228FastPrimitivePlace(plane,F,U,V,bc,values)
 BB,CC:=Explode(values);K:=Parent(BB);k:=PrimeField(K);p:=#k;assert CC ne 1;tt:=2*BB/(1-CC)-1;assert tt ne 0;qq:=((1-CC)*tt^2/2-1)/tt;R<W>:=PolynomialRing(K);candidates:=[];
 if qq ne 1 then uu:=W;vv:=((qq+1)*W+2)/(qq-1);else uu:=R!(-1);vv:=W;end if;
 nt:=uu^3+uu^2*vv+2*uu^2-uu*vv^2-vv^3-2*vv^2;dt:=uu^3+uu^2*vv+2*uu^2+uu*vv^2+2*uu*vv+vv^3+2*vv^2-2;ff:=nt-tt*dt;assert ff ne 0;
 for row in Roots(ff) do
  z:=row[1];uv:=[Evaluate(uu,z),Evaluate(vv,z)];ua,va:=Explode(uv);if ua eq va then continue;end if;
  if Evaluate(plane,[ua,va,K!1]) ne 0 then continue;end if;dv:=ua^3+ua^2*va+2*ua^2+ua*va^2+2*ua*va+va^3+2*va^2-2;if dv eq 0 then continue;end if;
  tq:=(ua^3+ua^2*va+2*ua^2-ua*va^2-va^3-2*va^2)/dv;qqq:=(-ua-va-2)/(ua-va);assert tq eq tt and qqq eq qq;assert (tt+1)*(qq*tt+1)/tt^2 eq BB and (tt^2-2*qq*tt-2)/tt^2 eq CC;Append(~candidates,uv);
 end for;assert #candidates eq 1;ua,va:=Explode(candidates[1]);eu:=Degree(MinimalPolynomial(ua,k));ev:=Degree(MinimalPolynomial(va,k));ee:=LCM(eu,ev);rr:=ee div eu;assert rr le 15;Y<Z>:=PolynomialRing(K);gg:=&*[Z-va^(p^(eu*j)):j in[0..rr-1]];M:=Matrix(k,[Eltseq(ua^i,k):i in[0..eu-1]]);assert Rank(M) eq eu;cc:=[];
 for z in Coefficients(gg) do assert z^(p^eu) eq z;sol:=Eltseq(Solution(M,Vector(k,Eltseq(z,k))));assert &+[K!sol[i+1]*ua^i:i in[0..eu-1]] eq z;Append(~cc,sol);end for;
 fu:=MinimalPolynomial(ua,k);gF:=&+[Evaluate(PolynomialRing(k)!cc[j+1],U)*V^j:j in[0..rr]];O:=MaximalOrderFinite(F);II:=ideal<O|Evaluate(fu,U),gF>;factors:=Factorization(II);ps:=[Place(z[1]):z in factors|z[2] gt 0];assert #ps eq 1 and Degree(ps[1]) eq ee;P:=ps[1];assert Valuation(Evaluate(fu,U),P) gt 0 and (gF eq 0 or Valuation(gF,P) gt 0);assert Y228PrimitiveJointKey([Evaluate(z,P):z in bc]) eq Y228PrimitiveJointKey(values);return P;
end function;
