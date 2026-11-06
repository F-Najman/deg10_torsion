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
a:=eval Read("work/deg10/continuation_20261029/M6N12_D5_all_source_lifts.m");lifts:=a[3];
X:=MDX612(GF(13));C:=Curve(X);FF:=FunctionField(C);AF,af:=AlgorithmicFunctionField(FF);B:=BaseField(AF);OF:=MaximalOrderFinite(AF);OI:=MaximalOrderInfinite(AF);
blocks,sigs,ds:=CuspBlocks(X);assert sigs eq a[1] and ds eq a[2];c0:=Support(blocks[Index(ds,1)])[1];cs:=Seqset(Cusps(X));
function PushPolynomial(D)
 result:=B!1;pp,mm:=Support(D);
 for i in [1..#pp] do
  I:=Ideal(pp[i]);
  if Order(I) eq OF then f:=Minimum(I);assert Degree(f) gt 0 and Degree(pp[i]) mod Degree(f) eq 0;result*:=(B!f)^(mm[i]*(Degree(pp[i]) div Degree(f)));
  else assert Order(I) eq OI;end if;
 end for;return result;
end function;
h0:=PushPolynomial(FunctionFieldDivisor(Divisor(c0)));hh:=[PushPolynomial(FunctionFieldDivisor(b)):b in blocks];
function Tuples(F,k)
 if k eq 0 then return [[]];end if;return [[a] cat v:a in F,v in Tuples(F,k-1)];
end function;
records:=[* *];left:=[];
for k in [k:k in [1..3360]|(k-1) mod 12 eq 5] do
 v:=lifts[k][2];L:=10*Divisor(c0)+&+[v[i]*(blocks[i]-ds[i]*Divisor(c0)):i in [1..24]];V,mp:=RiemannRochSpace(L);dim:=Dimension(V);rows:=[* *];positive:=false;
 if dim gt 4 then Append(~left,k);Append(~records,[*k,v,dim,rows,"dimension cap"*]);continue;end if;
 hl:=h0^10*&*[(hh[i]/h0^ds[i])^v[i]:i in [1..24]];
 for j in [1..dim] do for tail in Tuples(GF(13),dim-j) do
  u:=[GF(13)!0:i in [1..j-1]] cat [GF(13)!1] cat tail;f:=mp(V!u);z:=hl*(B!Norm(af(f)));assert Degree(Denominator(z)) eq 0;
  pol:=Numerator(z);pol/:=LeadingCoefficient(pol);assert Degree(pol) le 10;pass:=IsSquare(pol);degs:=[];mults:=[];even:=false;cusp:=false;
  if pass or #rows eq 0 then
   T:=L+Divisor(f);assert IsEffective(T) and Degree(T) eq 10;sp,mults:=Support(T);degs:=[Degree(P):P in sp];even:=&and[IsEven(m):m in mults];cusp:=&and[P in cs:P in sp];
   normT:=PushPolynomial(FunctionFieldDivisor(T));assert Degree(Denominator(normT)) eq 0;pp:=Numerator(normT);pp/:=LeadingCoefficient(pp);assert pp eq pol;
   if even then assert pass;end if;
  end if;
  Append(~rows,[*[Integers()!c:c in u],[Integers()!c:c in Coefficients(pol)],pass,degs,mults,even,cusp*]);
  if pass and even and not cusp then positive:=true;break;end if;
 end for;if positive then break;end if;end for;
 if positive then Append(~left,k);else assert #rows eq (13^dim-1) div 12;end if;
 Append(~records,[*k,v,dim,rows,positive*]);print "DECIC612_NORM_RR13_CLASS",k,dim,#rows,positive,Cputime(tt);
end for;
Write("work/deg10/continuation_20261029/cont_RR612_norm13_shard5_arithmetic.txt",Sprint([*records,left*]));print "DECIC612_NORM_RR13_COMPLETE",#records,left,Cputime(tt);quit;
