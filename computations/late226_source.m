SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(3*10^9);t0:=Cputime();
J:=JH(52,[27]);assert Dimension(J) eq 25;fp:=FrobeniusPolynomial(J,3);cs:=[1:i in [1..18]] cat [3:i in [1..6]];tr:=[];closed:=[];
for d in [1..10] do
 s:=-d*Coefficient(fp,50-d);for j in [1..d-1] do s-:=Coefficient(fp,50-j)*tr[d-j];end for;Append(~tr,Integers()!s);
 n:=3^d+1-s;for e in Divisors(d) do if e lt d then n-:=e*closed[e];end if;end for;assert n mod d eq 0;Append(~closed,n div d);
end for;
nc:=[closed[d]-#[e:e in cs|e eq d]:d in [1..10]];assert nc[1..4] eq [0,0,0,0];
Write("work/deg10/continuation_20261029b/C226/frobenius.m",Sprint([*Coefficients(fp),cs,nc*]));print "C226_FROBENIUS",nc,Cputime(t0);
for D in [5..10] do
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
 seen:={K|};orbs:={};points:=AssociativeArray();records:=[* *];jcount:=0;curves:=0;
 for j in K do
  if j in seen then continue;end if;for i in [0..D-1] do Include(~seen,j^(3^i));end for;jcount+:=1;
  es:=Twists(EllipticCurveWithjInvariant(j));
  for ee in es do
   curves+:=1;a:=aInvariants(ee);E:=EllipticCurve([K!0,a[2]+a[1]^2/4,K!0,a[4]+a[1]*a[3]/2,a[5]+a[3]^2/4]);
   if #E mod 52 ne 0 then continue;end if;
   ai:=aInvariants(E);cub:=x^3+ai[2]*x^2+ai[4]*x+ai[5];roots:=Roots(cub);if #roots ne 3 then continue;end if;
   ps:=[];subgroups:={};for rr in Roots(DivisionPolynomial(E,13)) do
    if not IsSquare(Evaluate(cub,rr[1])) then continue;end if;
    U:=E![rr[1],Sqrt(Evaluate(cub,rr[1]))];assert 13*U eq E!0;if U eq E!0 then continue;end if;
    ck:=Sort([Eltseq((k*U)[1]):k in [1..6]|GCD(k,13) eq 1]);if ck in subgroups then continue;end if;Include(~subgroups,ck);Append(~ps,U);
   end for;
   if #ps eq 0 then continue;end if;
   P:=E![roots[1][1],0];T:=E![roots[2][1],0];ts:=[P,T,P+T];
   for V in ps do
   ims:={};
   for a in [1..3],b in [1..3],k in [1..6] do if a eq b or GCD(k,13) ne 1 then continue;end if;Include(~ims,Key(BC(E,ts[a],ts[b]+k*V)));end for;
   ky:=Sort(Setseq(ims))[1];if ky in orbs then continue;end if;Include(~orbs,ky);
   vs:=[Evaluate(R!ky[1..D],K.1),Evaluate(R!ky[D+1..2*D],K.1)];d:=Deg(vs);
   assert &and[not IsDefined(points,z):z in ims];for z in ims do points[z]:=d;end for;Append(~records,[*d,ky*]);
   end for;
  end for;
 end for;
 counts:=[#[z:z in Keys(points)|points[z] eq d]:d in Divisors(D)];assert #seen eq #K and counts eq [nc[d]:d in Divisors(D)];
 Write(Sprintf("work/deg10/continuation_20261029b/C226/source%o.m",D),Sprint([*D,Coefficients(DefiningPolynomial(K,Fp)),true,jcount,curves,counts,records*]));
 print "C226_SOURCE_COMPLETE",D,counts,#records,Cputime(t0);
end for;
print "C226_SOURCE_ALL_COMPLETE",Cputime(t0);quit;
