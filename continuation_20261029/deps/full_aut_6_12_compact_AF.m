intrinsic Decic612OrbitsCompactAF(X::MDCrvMod11,xs::SeqEnum) -> SeqEnum,SeqEnum,SeqEnum,SeqEnum
{ Work entirely with native function-field places except for orbit representatives }
 p:=#BaseRing(X);d:=Degree(xs[1]);assert IsPrime(p);K:=GF(p^d);ix:=AssociativeArray();
 for i in [1..#xs] do assert not IsDefined(ix,xs[i]);ix[xs[i]]:=i;end for;
 C:=Curve(X);AF,af:=AlgorithmicFunctionField(FunctionField(C));O:=MaximalOrderFinite(AF);OI:=MaximalOrderInfinite(AF);
 function Recover(uv)
  f,g:=MDTwoGenerators(uv[1],uv[2]);ff:=af(Evaluate(f,X`_coordinates));gg:=af(Evaluate(g,X`_coordinates));I:=ideal<O|ff,gg>;fac:=Factorization(I);ps:=[Place(z[1]):z in fac|z[2] gt 0];
  if #ps eq 0 then fac:=Factorization(ideal<OI|ff,gg>);ps:=[Place(z[1]):z in fac|z[2] gt 0];end if;
  assert #ps eq 1 and Valuation(ff,ps[1]) gt 0 and Valuation(gg,ps[1]) gt 0;return ps[1];
 end function;
 keys:=[[]:x in xs];labels:=[0:x in xs];rs:=[];repi:=[];tt:=Cputime();
 for i in [1..#xs] do
  if labels[i] ne 0 then continue;end if;
  x:=CurvePlace(C,xs[i]);Append(~rs,x);Append(~repi,i);ind:=#rs;E0:=EllipticCurve(X,x);ls:=LevelStructure612(X,x);
  vals:=MDM2ToCommonField(aInvariants(E0) cat [ls`P[1],ls`P[2],ls`Q[1],ls`Q[2]],K);E:=EllipticCurve(vals[1..5]);P:=E!vals[6..7];Q:=E!vals[8..9];assert Order(P) eq 6 and Order(Q) eq 12;
  assert WeilPairing(P,2*Q,6) eq K!(PrimitiveElement(BaseRing(X))^((p-1) div 6));
  for par in Decic612Parameters() do
   a,b,c,dd:=Explode(par);PP:=a*P+2*b*Q;QQ:=c*P+dd*Q;uv:=Decic612Normalize(E,PP,QQ);im:=Recover(uv);assert Degree(im) eq d and IsDefined(ix,im);j:=ix[im];ky:=MDM2FrobeniusKey(uv,K,p,d);
   if labels[j] eq 0 then labels[j]:=ind;keys[j]:=ky;else assert labels[j] eq ind and keys[j] eq ky;end if;
  end for;assert labels[i] eq ind;
  if ind mod 25 eq 0 then print "DECIC612_COMPACT_AF_AUT_PROGRESS",ind,Cputime(tt);end if;
 end for;
 assert &and[a gt 0:a in labels] and #Seqset(keys) eq #xs;return rs,labels,keys,repi;
end intrinsic;
