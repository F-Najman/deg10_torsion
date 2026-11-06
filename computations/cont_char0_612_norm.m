SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(8*10^9);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
data:=eval Read("work/deg10/continuation_20261029/M6N12_D5_all_source_lifts.m");
R<T>:=PolynomialRing(Rationals());F<z>:=NumberField(T^2-T+1);
 A<u,v>:=AffineSpace(F,2);
 rawmodel:=Split(Read("work/deg10/shared/mdmagma/models_X1_m_n/X1_6_12.txt"));
 for line in rawmodel do
  val:=Split(Split(line,"=")[2],";")[1];
  if line[1] eq "X" then equation:=eval(val);end if;
  if line[1] eq "r" then r:=eval(val);end if;
  if line[1] eq "s" then s:=eval(val);end if;
  if line[1] eq "q" then q:=eval(val);end if;
  if line[1] eq "E" then es:=eval(val);end if;
 end for;
 C:=ProjectiveClosure(Curve(A,equation));FF:=FunctionField(C);
 X:=New(MDCrvMod11);X`M:=6;X`N:=12;
 X`_P:=[FF!q];X`_Q:=[FF!0,FF!0];X`_coordinates:=[FF.1,FF.2];
 _InitMDCrvMod(X,12,C,[FF!a:a in es]);

C:=Curve(X);FF:=FunctionField(C);print "M6N12_CHAR0_MODEL",Genus(X),Cputime(tt);
 E:=EllipticCurve([FF!a:a in X`_E]);Q:=E![0,0];assert 12*Q eq E!0;
 px:=FF!X`_P[1];ai:=aInvariants(E);R<Y>:=PolynomialRing(FF);
 ys:=Roots(Y^2+(ai[1]*px+ai[3])*Y-px^3-ai[2]*px^2-ai[4]*px-ai[5]);
 assert #ys eq 2;P:=E![px,ys[1][1]];

 if WeilPairing(P,2*Q,6) ne z then P:=-P;end if;
 assert WeilPairing(P,2*Q,6) eq z;
 print "M6N12_CHAR0_WEIL",Cputime(tt);
 fs:=[(P+5*Q)[1]-P[1],(2*P+5*Q)[1]-P[1],(P+2*Q)[1]];
 assert &and[f ne 0:f in fs];print "M6N12_CHAR0_THREE_UNITS",Cputime(tt);
 cs:=Cusps(X);groups:=AssociativeArray();
 for c in cs do
  s:=[Valuation(f,c):f in fs];
  if IsDefined(groups,s) then Append(~groups[s],c);else groups[s]:=[c];end if;
 end for;
 cols:=[34,58,27];profiles:=[[sig[j]:j in cols]:sig in data[1]];
 assert #Seqset(profiles) eq 24 and Seqset(profiles) eq Keys(groups);
 blocks:=[&+[Divisor(c):c in groups[sig]]:sig in profiles];ds:=[Degree(D):D in blocks];
 print "M6N12_CHAR0_CUSP_BLOCKS",#blocks,ds,Cputime(tt);
 sigs:=data[1];
 c0:=Support(blocks[Index(ds,1)])[1];

index:=2791;effective:=eval Read("work/deg10/continuation_20261029/M6N12_effective_quintic_lifts_v2.m");coeff:=effective[1][2];L:=&+[coeff[i]*blocks[i]:i in [1..24]];assert Degree(L) eq 10 and IsEffective(L);
V,mp:=RiemannRochSpace(L);assert Dimension(V) eq 2;print "DECIC612_CHAR0_NORM_RR_READY",Cputime(tt);
for j in [1..2] do
 f:=mp(V.j);vs:=[Valuation(f,c):c in cs];if &or[e lt 0:e in vs] then break;end if;
end for;
poles:=&+[-vs[i]*Divisor(cs[i]):i in [1..#cs]|vs[i] lt 0];degree:=Degree(poles);assert degree le 10;
print "DECIC612_CHAR0_NORM_POLE_DEGREE",degree,Cputime(tt);
assert degree eq 10 and poles eq L;
AF,af:=AlgorithmicFunctionField(FF);B:=BaseField(AF);OF:=MaximalOrderFinite(AF);OI:=MaximalOrderInfinite(AF);IL,JL:=Ideals(FunctionFieldDivisor(L));assert IsIntegral(IL*(af(f)*OF)) and IsIntegral(JL*(af(f)*OI));
function PushPolynomial(D)
 result:=B!1;pp,mm:=Support(D);
 for i in [1..#pp] do I:=Ideal(pp[i]);if Order(I) eq OF then h:=Minimum(I);assert Degree(h) gt 0 and Degree(pp[i]) mod Degree(h) eq 0;result*:=(B!h)^(mm[i]*(Degree(pp[i]) div Degree(h)));else assert Order(I) eq OI;end if;end for;return result;
end function;
hL:=PushPolynomial(FunctionFieldDivisor(L));gen:=AF.1;n:=Degree(AF);assert n eq 6;
rows:=[];
for i in [0..n-1] do
 product:=af(f)*gen^i;vv:=Eltseq(product);assert #vv eq n and &+[(AF!vv[j])*gen^(j-1):j in [1..n]] eq product;Append(~rows,vv);
end for;
M:=Matrix(B,rows);cp:=CharacteristicPolynomial(M);assert Degree(cp) eq 6;print "DECIC612_CHAR0_NORM_POLYNOMIAL_READY",Cputime(tt);
Ra<a>:=PolynomialRing(F);hu:=[Ra!0:j in [0..10]];
for k in [0..6] do
 rat:=hL*(B!Coefficient(cp,k));assert Degree(Denominator(rat)) eq 0;poly:=Numerator(rat)/Coefficient(Denominator(rat),0);assert Degree(poly) le 10;
 for j in [0..Degree(poly)] do hu[j+1]+:=Coefficient(poly,j)*a^k;end for;
end for;lead:=hu[11];assert lead ne 0;Rat:=FieldOfFractions(Ra);cc:=[Rat!q/lead:q in hu];bb:=[Rat!0:j in [0..5]];bb[6]:=1;
for j in [4..0 by -1] do total:=Rat!0;for i in [j+1..4] do total+:=bb[i+1]*bb[5+j-i+1];end for;bb[j+1]:=(cc[5+j+1]-total)/2;end for;
G:=Ra!0;residuals:=[];
for j in [0..4] do
 total:=Rat!0;for i in [0..j] do total+:=bb[i+1]*bb[j-i+1];end for;rr:=cc[j+1]-total;Append(~residuals,rr);G:=GCD(G,Numerator(rr));
end for;assert G ne 0;G/:=LeadingCoefficient(G);candidate:=G*lead;candidate/:=LeadingCoefficient(candidate);
roots:=Roots(candidate);print "DECIC612_CHAR0_NORM_BRANCH_GCD",G,"LEAD",lead,"F_ROOTS",roots,Cputime(tt);
function EncF(x)
 cc:=Eltseq(F!x);assert #cc eq 2;return [[Numerator(q),Denominator(q)]:q in cc];
end function;
function EncP(p)
 return [EncF(x):x in Coefficients(p)];
end function;
function EncR(r)
 return [*EncP(Numerator(B!r)),EncP(Denominator(B!r))*];
end function;
fibers:=[* *];found:=false;
for root in roots do
 value:=root[1];Rpol<uu>:=PolynomialRing(F);pol:=Rpol![Evaluate(h,value):h in hu];assert pol ne 0;pol/:=LeadingCoefficient(pol);pass:=IsSquare(pol);
 degs:=[];mult:=[];
 if pass then
  T:=poles+Divisor(f-value);assert IsEffective(T) and Degree(T) eq 10;sp,mult:=Support(T);degs:=[Degree(P):P in sp];
  if #sp eq 1 and degs eq [5] and mult eq [2] then
   P:=sp[1];Ej:=EllipticCurve(X,P);found:=true;Write("work/deg10/continuation_20261029/M6N12_quintic_occurrence_candidate.txt",Sprint([*P,aInvariants(Ej),jInvariant(Ej)*]));print "DECIC612_CHAR0_NORM_QUINTIC_FOUND",value;
  end if;
 end if;
 Append(~fibers,[*EncF(value),EncP(pol),pass,degs,mult*]);
end for;
out:=[*index,degree,coeff,[EncR(x):x in Coefficients(DefiningPolynomial(AF))],[EncR(x):x in Eltseq(af(f))],EncR(hL),[EncP(h):h in hu],EncP(G),EncP(lead),fibers,found*];
Write("work/deg10/continuation_20261029/M6N12_char0_norm_certificate.m",Sprint(out));Write("work/deg10/continuation_20261029/M6N12_char0_norm_function.m",Sprint(f));print "DECIC612_CHAR0_NORM_COMPLETE",found,Cputime(tt);quit;
