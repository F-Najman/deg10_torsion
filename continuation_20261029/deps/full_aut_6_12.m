intrinsic Decic612Parameters() -> SeqEnum
{ Pairing-preserving changes of the (6,12) level structure, modulo simultaneous sign }
 pars:=[];
 for a in [1,3,5] do for b,c in [0..5] do for d in [1,3,5,7,9,11] do
  if (a*d-2*b*c) mod 6 ne 1 then continue;end if;
  par:=[a,b,c,d];neg:=[(-a) mod 6,(-b) mod 6,(-c) mod 6,(-d) mod 12];
  if par lt neg then Append(~pars,par);end if;
 end for;end for;end for;assert #pars eq 96;return pars;
end intrinsic;
intrinsic Decic612Normalize(E::CrvEll,P::PtEll,Q::PtEll) -> SeqEnum
{ Exact u,v inverse from the accepted model_6_12 normalization }
 assert Order(P) eq 6 and Order(Q) eq 12 and Order(WeilPairing(P,2*Q,6)) eq 6;
 a1,a2,a3,a4,a6:=Explode(aInvariants(E));x0:=Q[1];y0:=Q[2];
 aa3:=2*y0+a3+a1*x0;aa2:=3*x0+a2;aa4:=3*x0^2+2*x0*a2+a4-a1*y0;slope:=aa4/aa3;
 aaa1:=2*slope+a1;aaa2:=aa2-slope^2-a1*slope;scale:=aaa2/aa3;
 b:=-aaa2^3/aa3^2;c:=1-aaa1*scale;r:=b/c;s:=c^2/(b-c);u:=(r-1)/(r-s);
 q:=scale^2*(P[1]-x0);v:=q*u^2*(u-1)^2/(u+1)-u^2-1;
 assert r eq (u^2+1)/(u^2-u) and s eq (u^2+u+1)/u^2;return [u,v];
end intrinsic;
intrinsic Decic612Orbits(X::MDCrvMod11,xs::SeqEnum) -> SeqEnum,SeqEnum,SeqEnum
{ Exhaustive verified finite-field orbit partition under the 96 maps }
 p:=#BaseRing(X);d:=Degree(xs[1]);K:=GF(p^d);keys:=[MDM2FrobeniusKey([Evaluate(f,x):f in X`_coordinates],K,p,d):x in xs];
 allowed:=Seqset(keys);assert #allowed eq #xs;seen:=AssociativeArray();rs:=[];labels:=[];
 for i in [1..#xs] do
  ky:=keys[i];
  if not IsDefined(seen,ky) then
   x:=xs[i];Append(~rs,x);j:=#rs;E:=EllipticCurve(X,x);ls:=LevelStructure612(X,x);P:=ls`P;Q:=ls`Q;
   for par in Decic612Parameters() do
    a,b,c,dd:=Explode(par);PP:=a*P+2*b*Q;QQ:=c*P+dd*Q;
    im:=MDM2FrobeniusKey(Decic612Normalize(E,PP,QQ),K,p,d);assert im in allowed;
    if IsDefined(seen,im) then assert seen[im] eq j;else seen[im]:=j;end if;
   end for;
  end if;
  Append(~labels,seen[ky]);
 end for;assert #Keys(seen) eq #xs;return rs,labels,keys;
end intrinsic;
