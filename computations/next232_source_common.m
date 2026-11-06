SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(3*10^9);t0:=Cputime();
data:=eval Read("work/deg10/C232_resolution_20261030/frobenius.m");nc:=data[3];
for D in DS do
 K:=GF(3^D);Fp:=GF(3);R<x>:=PolynomialRing(K);
 function BC(E,P,Q)
  ai:=aInvariants(E);xx:=Q[1]-P[1];yy:=Q[2];return [xx^3/yy^2,(ai[2]+3*P[1])*xx^2/yy^2];
 end function;
 function Key(v)
  keys:=[];for i in [1..D] do Append(~keys,&cat[[Integers()!a:a in Eltseq(z,Fp)]:z in v]);v:=[z^3:z in v];end for;return Sort(keys)[1];
 end function;
 function Deg(v)
  for d in Divisors(D) do if &and[z^(3^d) eq z:z in v] then return d;end if;end for;error "fielddegree";
 end function;
 seen:={K|};orbs:={};points:=AssociativeArray();records:=[* *];jcount:=0;curves:=0;torsioncurves:=0;
 for j in K do
  if j in seen then continue;end if;for i in [0..D-1] do Include(~seen,j^(3^i));end for;jcount+:=1;
  for ee in Twists(EllipticCurveWithjInvariant(j)) do
   curves+:=1;a:=aInvariants(ee);E:=EllipticCurve([K!0,a[2]+a[1]^2/4,K!0,a[4]+a[1]*a[3]/2,a[5]+a[3]^2/4]);
   if #E mod 64 ne 0 then continue;end if;
   ai:=aInvariants(E);cub:=x^3+ai[2]*x^2+ai[4]*x+ai[5];if #Roots(cub) ne 3 then continue;end if;
   G,gm:=AbelianGroup(E);S:=sub<G|[(Order(G.i) div 2^Valuation(Order(G.i),2))*G.i:i in [1..Ngens(G)]]>;
   qs:=[G!q:q in S|Order(q) eq 32];if #qs eq 0 then continue;end if;ps:=[G!p:p in S|Order(p) eq 2];torsioncurves+:=1;hseen:={};
   for q in qs,p in ps do
    if p eq 16*q then continue;end if;H:=sub<G|q,p>;assert #H eq 64;hk:=&cat Sort([Eltseq(G!h):h in H]);if hk in hseen then continue;end if;Include(~hseen,hk);
    P:=gm(p);Q:=gm(q);assert Order(P) eq 2 and Order(Q) eq 32 and P ne 16*Q;ims:={};
    for b,c in [0,1],a in [1,3,5,7,9,11,13,15] do Include(~ims,Key(BC(E,P+b*(16*Q),c*P+a*Q)));end for;
    ky:=Sort(Setseq(ims))[1];if ky in orbs then continue;end if;Include(~orbs,ky);
    vs:=[Evaluate(R!ky[1..D],K.1),Evaluate(R!ky[D+1..2*D],K.1)];d:=Deg(vs);
    assert &and[not IsDefined(points,z):z in ims];for z in ims do points[z]:=d;end for;Append(~records,[*d,ky*]);
   end for;
  end for;
  if jcount mod 500 eq 0 then print "C232_SOURCE_PROGRESS",D,jcount,#records,Cputime(t0);end if;
 end for;
 counts:=[#[z:z in Keys(points)|points[z] eq d]:d in Divisors(D)];print "C232_SOURCE_COUNTS",D,counts,[nc[d]:d in Divisors(D)];assert #seen eq #K and counts eq [nc[d]:d in Divisors(D)];
 Write(Sprintf("work/deg10/C232_resolution_20261030/source%o.m",D),Sprint([*D,Coefficients(DefiningPolynomial(K,Fp)),true,jcount,curves,counts,records*]):Overwrite:=true);
 print "C232_SOURCE_COMPLETE",D,counts,#records,torsioncurves,Cputime(t0);
end for;
print "C232_SOURCE_BATCH_COMPLETE",Cputime(t0);quit;
