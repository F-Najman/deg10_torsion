SetColumns(0);SetNthreads(1);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
function CuspUnitsM4(X)
 FF:=FunctionField(Curve(X));E:=EllipticCurve([FF!f:f in X`_E]);P:=E![FF!f:f in X`_P];Q:=E![FF!f:f in X`_Q];
 ai:=aInvariants(E);B:=ai[2];assert ai eq [FF!1,B,B,FF!0,FF!0] and P eq E![0,0];
 assert Discriminant(E) eq B^4*(1-16*B);fs:=[B,1-16*B];
 for k in [1..8] do for j in [0..3] do
  A:=k*Q+j*P;assert A ne E!0 and A[1] ne 0 and A[2] ne 0 and A[1]+B ne 0;
  fs cat:=[A[1],A[2],A[1]+B,A[2]+A[1]+B];
 end for;end for;
 for k in [1..8] do A:=k*Q;for j in [1..3] do R:=A+j*P;assert A ne R and A ne -R;Append(~fs,A[1]-R[1]);end for;end for;
 rawcs:=Cusps(X);groups:=AssociativeArray();
 for c in rawcs do sig:=[Valuation(f,c):f in fs];
  if IsDefined(groups,sig) then Append(~groups[sig],c);else groups[sig]:=[c];end if;
 end for;
 sigs:=Sort(Setseq(Keys(groups)));assert #sigs eq 24;
 return [&+[Divisor(c):c in groups[sig]]:sig in sigs],sigs;
end function;
a:=eval Read("work/deg10/continuation_20261029/diagonal_4_16_D5_arithmetic.txt");
FB:=FreeAbelianGroup(#a[7][2][1]);G,pi:=quo<FB|[FB!r:r in a[6]]>;
Free:=FreeAbelianGroup(24);h:=hom<Free->G|[pi(FB!v):v in a[7][2]]>;
assert Invariants(Image(h)) eq [2,2,2,4,4,20,40,40,80];
KB:=Matrix(Integers(),[Eltseq(Free!w):w in Generators(Kernel(h))]);LL:=LatticeWithBasis(LLL(KB));
classes:=AssociativeArray();lifts:=[* *];
for mm in a[10] do
 target:=a[9][mm[1]];cl:=2*pi(FB!target[4]);key:=Eltseq(cl);
 if IsDefined(classes,key) then Append(~lifts[classes[key]][1],mm[1]);continue;end if;
 vv:=Eltseq(cl@@h);w:=ClosestVector(LL,Vector(Rationals(),vv));v:=[vv[i]-Integers()!w[i]:i in [1..24]];assert h(Free!v) eq cl;
 Append(~lifts,[*[mm[1]],v,key*]);classes[key]:=#lifts;
end for;
print "DECIC416_RR13_CLASS_COUNT",#lifts,Cputime(t0);
Write("work/deg10/continuation_20261029/M4N16_D5_cusp_lifts.m",Sprint([*a[1],a[2],lifts*]));
X:=MDX11(4,16,GF(13));C:=Curve(X);blocks,sigs:=CuspUnitsM4(X);ds:=[Degree(b):b in blocks];assert sigs eq a[1] and ds eq a[2];c0:=Support(blocks[Index(ds,1)])[1];
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
 for j in [1..dim] do for tail in Tuples(GF(13),dim-j) do
  u:=[GF(13)!0:i in [1..j-1]] cat [GF(13)!1] cat tail;f:=mp(V!u);T:=L+Divisor(f);assert IsEffective(T) and Degree(T) eq 10;sp,ms:=Support(T);degrees:=[Degree(P):P in sp];even:=&and[IsEven(m):m in ms];
  Append(~rows,[*[Integers()!c:c in u],degrees,ms,even*]);
  if even then positive:=true;break;end if;
 end for;if positive then break;end if;end for;
 if positive then Append(~remaining,k);else assert #rows eq (13^dim-1) div 12;end if;
 Append(~records,[*k,lifts[k][1],v,dim,rows,positive*]);
 print "DECIC416_RR13_PARITY",k,#lifts,"dimension",dim,"tested",#rows,"even_member",positive,Cputime(t0);
 Write("work/deg10/continuation_20261029/M4N16_D5_RR13_partial.m",Sprint([*records,remaining*]):Overwrite:=true);
end for;
Write("work/deg10/continuation_20261029/M4N16_D5_RR13_arithmetic.txt",Sprint([*records,remaining*]));
print "DECIC416_RR13_PARITY_COMPLETE",#lifts,remaining,Cputime(t0);quit;
