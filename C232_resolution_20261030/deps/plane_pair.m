/* Invert the supplied q,t plane parametrization before constructing a place.
Every candidate is checked on the exact plane model and both moduli coordinates. */
function C232PlanePair(C,bc)
 K:=Parent(bc[1]);b,c:=Explode(bc);assert b ne 0 and c ne 1;t:=2*b/(1-c)-1;assert t ne 0;q:=((1-c)*t^2-2)/(2*t);
 R<v>:=PolynomialRing(K);L:=FieldOfFractions(R);vv:=L!v;uu:=vv*(q+1)/(q-1-2*vv);
 nn:=uu^3*vv+uu^3/2+uu^2*vv/2-uu*vv^3-uu*vv^2/2-vv^3/2;
 dd:=uu^3*vv^3-uu^3*vv-uu^3/2-uu^2*vv^2-uu^2*vv/2-uu*vv^3-uu*vv^2/2-vv^3/2;
 pol:=R!Numerator(t*dd-nn);assert pol ne 0;solutions:=[];
 for rr in Roots(pol) do
  v0:=rr[1];den:=q-1-2*v0;if den eq 0 then continue;end if;u0:=v0*(q+1)/den;if u0 eq v0 then continue;end if;
  nd:=u0^3*v0+u0^3/2+u0^2*v0/2-u0*v0^3-u0*v0^2/2-v0^3/2;
  de:=u0^3*v0^3-u0^3*v0-u0^3/2-u0^2*v0^2-u0^2*v0/2-u0*v0^3-u0*v0^2/2-v0^3/2;
  if de eq 0 or nd/de ne t then continue;end if;
  qq:=(2*u0*v0+u0+v0)/(u0-v0);if qq ne q then continue;end if;
  if Evaluate(DefiningPolynomial(C),[u0,v0,K!1]) ne 0 then continue;end if;
  assert (t+1)*(qq*t+1)/t^2 eq b and (t^2-2*qq*t-2)/t^2 eq c;Append(~solutions,[u0,v0]);
 end for;
 assert #solutions eq 1;return solutions[1];
end function;
