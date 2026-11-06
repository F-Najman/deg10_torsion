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
a:=eval Read("work/deg10/continuation_20261029/M6N12_D5_all_source_lifts.m");lifts:=a[3];assert #lifts eq 3360;
X:=MDX612(GF(13));C:=Curve(X);blocks,sigs,ds:=CuspBlocks(X);assert sigs eq a[1] and ds eq a[2];c0:=Support(blocks[Index(ds,1)])[1];cs:=Seqset(Cusps(X));
function Tuples(F,k)
 if k eq 0 then return [[]];end if;return [[a] cat v:a in F,v in Tuples(F,k-1)];
end function;
records:=[* *];left:=[];
for k in [1..#lifts] do
 if (k-1) mod 12 ne 2 then continue;end if;
 v:=lifts[k][2];L:=10*Divisor(c0)+&+[v[i]*(blocks[i]-ds[i]*Divisor(c0)):i in [1..24]];assert Degree(L) eq 10;
 V,mp:=RiemannRochSpace(L);dim:=Dimension(V);assert dim ge 2;rows:=[* *];positive:=false;
 if dim gt 3 then Append(~left,k);Append(~records,[*k,v,dim,rows,"dimension cap"*]);continue;end if;
 for j in [1..dim] do for tail in Tuples(GF(13),dim-j) do
  u:=[GF(13)!0:i in [1..j-1]] cat [GF(13)!1] cat tail;T:=L+Divisor(mp(V!u));assert IsEffective(T) and Degree(T) eq 10;sp,mults:=Support(T);degrees:=[Degree(P):P in sp];even:=&and[IsEven(m):m in mults];cusp:=&and[P in cs:P in sp];
  Append(~rows,[*[Integers()!c:c in u],degrees,mults,even,cusp*]);if even and not cusp then positive:=true;break;end if;
 end for;if positive then break;end if;end for;
 if positive then Append(~left,k);else assert #rows eq (13^dim-1) div 12;end if;
 Append(~records,[*k,v,dim,rows,positive*]);print "DECIC612_SOURCE_RR13_CLASS",k,dim,#rows,positive,Cputime(tt);
end for;
assert #records eq 280;Write("work/deg10/continuation_20261029/cont_RR612_source_p13_shard2_arithmetic.txt",Sprint([*records,left*]));print "DECIC612_SOURCE_RR13_SHARD_COMPLETE",2,#left,left,Cputime(tt);quit;
