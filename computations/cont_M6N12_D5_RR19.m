SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(8*10^9);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
Attach("work/deg10/continuation_20261029/deps/fast_hecke.m");
Attach("work/deg10/continuation_20261029/deps/fast_hecke_11.m");
Attach("work/deg10/continuation_20261029/deps/model_6_12.m");
function CuspBlocks(X)
C:=Curve(X);FF:=FunctionField(C);E:=EllipticCurve(X`_E);Q:=E![0,0];px:=FF!X`_P[1];ai:=aInvariants(E);R<Y>:=PolynomialRing(FF);
ys:=Roots(Y^2+(ai[1]*px+ai[3])*Y-px^3-ai[2]*px^2-ai[4]*px-ai[5]);assert #ys eq 2;P:=E![px,ys[1][1]];
z:=PrimitiveElement(BaseRing(X))^((#BaseRing(X)-1) div 6);if WeilPairing(P,2*Q,6) ne z then P:=-P;end if;assert WeilPairing(P,2*Q,6) eq z;
fs:=[-ai[2],1-ai[1],Discriminant(E)];anchor:=P[1];
for i in [0..5] do for j in [0..11] do
 A:=i*P+j*Q;if A eq E!0 then continue;end if;
 if A ne Q and A ne -Q then assert A[1] ne 0;Append(~fs,A[1]);end if;
 if A ne P and A ne -P then assert A[1] ne anchor;Append(~fs,A[1]-anchor);end if;
end for;end for;
cs:=Cusps(X);groups:=AssociativeArray();
for c in cs do
 sig:=[Valuation(f,c):f in fs];if IsDefined(groups,sig) then Append(~groups[sig],c);else groups[sig]:=[c];end if;
end for;
sigs:=Sort(Setseq(Keys(groups)));blocks:=[&+[Divisor(c):c in groups[sig]]:sig in sigs];ds:=[Degree(b):b in blocks];assert #blocks eq 24;
return blocks,sigs,ds;
end function;
a:=eval Read("work/deg10/continuation_20261029/diagonal_6_12_D5_arithmetic.txt");
FB:=FreeAbelianGroup(#a[5][2][1]);G,pi:=quo<FB|[FB!r:r in a[4][2]]>;
Free:=FreeAbelianGroup(24);h:=hom<Free->G|[pi(FB!v):v in a[5][2]]>;
assert #Image(h) eq 5971968;
KB:=Matrix(Integers(),[Eltseq(Free!w):w in Generators(Kernel(h))]);LL:=LatticeWithBasis(LLL(KB));
classes:=AssociativeArray();lifts:=[* *];
for mm in a[13] do
 target:=a[12][mm[1]];cl:=2*pi(FB!target[2]);key:=Eltseq(cl);
 if IsDefined(classes,key) then Append(~lifts[classes[key]][1],mm[1]);continue;end if;
 vv:=Eltseq(cl@@h);w:=ClosestVector(LL,Vector(Rationals(),vv));v:=[vv[i]-Integers()!w[i]:i in [1..24]];assert h(Free!v) eq cl;
 Append(~lifts,[*[mm[1]],v,key*]);classes[key]:=#lifts;
end for;
print "DECIC612_RR19_CLASS_COUNT",#lifts,Cputime(tt);
Write("work/deg10/continuation_20261029/M6N12_D5_cusp_lifts.m",Sprint([*a[1],a[2],lifts*]));
X:=MDX612(GF(19));C:=Curve(X);blocks,sigs,ds:=CuspBlocks(X);assert sigs eq a[1] and ds eq a[2];c0:=Support(blocks[Index(ds,1)])[1];
function Tuples(F,k)
 if k eq 0 then return [[]];end if;return [[a] cat v:a in F,v in Tuples(F,k-1)];
end function;
records:=[* *];remaining:=[];
for k in [1..#lifts] do
 v:=lifts[k][2];L:=10*Divisor(c0)+&+[v[i]*(blocks[i]-ds[i]*Divisor(c0)):i in [1..24]];assert Degree(L) eq 10;
 V,mp:=RiemannRochSpace(L);dim:=Dimension(V);rows:=[* *];positive:=false;
 if dim gt 3 then
  Append(~remaining,k);Append(~records,[*k,lifts[k][1],v,dim,"dimension exceeds enumeration cap"*]);continue;
 end if;
 for j in [1..dim] do for tail in Tuples(GF(19),dim-j) do
  u:=[GF(19)!0:i in [1..j-1]] cat [GF(19)!1] cat tail;f:=mp(V!u);T:=L+Divisor(f);assert IsEffective(T) and Degree(T) eq 10;sp,ms:=Support(T);degrees:=[Degree(P):P in sp];even:=&and[IsEven(m):m in ms];
  Append(~rows,[*[Integers()!c:c in u],degrees,ms,even*]);
  if even then positive:=true;break;end if;
 end for;if positive then break;end if;end for;
 if positive then Append(~remaining,k);else assert #rows eq (19^dim-1) div 18;end if;
 Append(~records,[*k,lifts[k][1],v,dim,rows,positive*]);
 print "DECIC612_RR19_PARITY",k,#lifts,"dimension",dim,"tested",#rows,"even_member",positive,Cputime(tt);
 Write("work/deg10/continuation_20261029/M6N12_D5_RR19_partial.m",Sprint([*records,remaining*]):Overwrite:=true);
end for;
Write("work/deg10/continuation_20261029/M6N12_D5_RR19_arithmetic.txt",Sprint([*records,remaining*]));
print "DECIC612_RR19_PARITY_COMPLETE",#lifts,remaining,Cputime(tt);quit;
