intrinsic Decic420Normalize(E::CrvEll,P::PtEll,Q::PtEll) -> SeqEnum
{ Canonical q,t coordinates for an independent level(4,20) pair }
 a1,a2,a3,a4,a6:=Explode(aInvariants(E));r:=P[1];t0:=P[2];
 s:=(3*r^2+2*a2*r+a4-a1*t0)/(2*t0+a1*r+a3);
 a1p:=a1+2*s;a2p:=a2-s*a1+3*r-s^2;a3p:=a3+r*a1+2*t0;
 a4p:=a4-s*a3+2*r*a2-(t0+r*s)*a1+3*r^2-2*s*t0;
 a6p:=a6+r*a4+r^2*a2+r^3-t0*a3-r*t0*a1-t0^2;
 assert a4p eq 0 and a6p eq 0 and a3p eq a1p*a2p and a1p ne 0;
 scale:=1/a1p;b:=scale^2*a2p;
 xx:=scale^2*(Q[1]-r);yy:=scale^3*(Q[2]-t0-s*(Q[1]-r));
 q:=1+2*b/xx;t:=1+4*yy/((q+1)*xx);
 assert b eq (q^2-1)*(t^2-1)/16 and xx eq (q+1)*(t^2-1)/8;
 assert yy eq (q+1)^2*(t-1)^2*(t+1)/32;
 return [q,t];
end intrinsic;
intrinsic Decic420Parameters() -> SeqEnum
{ All 96 symplectic changes of a (4,20) level pair, modulo simultaneous sign }
 ps:=[];
 for a,b,c in [0..3] do for d in [0..19] do
  if GCD(d,5) ne 1 or (a*d-5*b*c) mod 4 ne 1 then continue;end if;
  v:=[a,b,c,d];w:=[(-a) mod 4,(-b) mod 4,(-c) mod 4,(-d) mod 20];
  if v lt w then Append(~ps,v);end if;
 end for;end for;assert #ps eq 96;return ps;
end intrinsic;
intrinsic Decic420Orbits(X::MDCrvMod11,xs::SeqEnum) -> SeqEnum,SeqEnum,SeqEnum
{ Exhaustive actual place-key partition, fixing the quadratic residue base }
 F:=BaseRing(X);dd:=Degree(xs[1]);K:=ext<F|dd>;
 keys:=[M4SixteenRelativeKey([Evaluate(f,x):f in X`_coordinates],K,F,dd):x in xs];
 allowed:=Seqset(keys);assert #allowed eq #xs;seen:=AssociativeArray();rs:=[];labels:=[];
 for i in [1..#xs] do
  ky:=keys[i];
  if not IsDefined(seen,ky) then
   Append(~rs,xs[i]);j:=#rs;
   q,t:=Explode(M4SixteenRelativeToCommon([Evaluate(f,xs[i]):f in X`_coordinates],K,F));
   B:=(q^2-1)*(t^2-1)/16;E:=EllipticCurve([K!1,B,B,K!0,K!0]);
   P:=E![0,0];Q:=E![(q+1)*(t^2-1)/8,(q+1)^2*(t-1)^2*(t+1)/32];
   assert Order(P) eq 4 and Order(Q) eq 20 and Order(WeilPairing(P,5*Q,4)) eq 4;
   for par in Decic420Parameters() do
    a,b,c,d:=Explode(par);PP:=a*P+5*b*Q;QQ:=c*P+d*Q;
    assert Order(PP) eq 4 and Order(QQ) eq 20 and WeilPairing(PP,5*QQ,4) eq WeilPairing(P,5*Q,4);
    im:=M4SixteenRelativeKey(Decic420Normalize(E,PP,QQ),K,F,dd);assert im in allowed;
    if IsDefined(seen,im) then assert seen[im] eq j;else seen[im]:=j;end if;
   end for;
  end if;
  Append(~labels,seen[ky]);
 end for;assert #Keys(seen) eq #xs;return rs,labels,keys;
end intrinsic;
