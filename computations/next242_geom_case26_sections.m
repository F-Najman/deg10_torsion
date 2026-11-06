SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(4*10^9);tt:=Cputime();AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
model:=eval Read("work/deg10/next_cases_20261030/C242/geometry/diagonal17_normalized_model_data.m");param:=eval Read("work/deg10/next_cases_20261030/C242/geometry/diagonal17_parameter_data.m");
function MakeField(k,data)
 A<z>:=FunctionField(k);P:=PolynomialRing(k);R<X>:=PolynomialRing(A);function Rat(v)return A!(P!v[1])/A!(P!v[2]);end function;f:=R![Rat(v):v in data[1]];f/:=LeadingCoefficient(f);F<r>:=FunctionField(f);coords:=[F![Rat(v):v in row]:row in data[2]];return F,coords;
end function;
k5:=GF(5);F5,c5:=MakeField(k5,model);assert Genus(F5) eq 17 and Degree(F5) eq 8;P5:=PolynomialRing(k5);u5:=Evaluate(PolynomialRing(F5)!param[2][1],c5[1])/Evaluate(PolynomialRing(F5)!param[2][2],c5[1]);t5:=u5^2;D5:=t5^2+13*t5+49;N5:=t5^2+5*t5+1;j5:=N5^3*D5/t5;cs5:=Poles(j5);assert Sort([Degree(p):p in cs5]) eq [1,1,1,1,1,1,2,2,2,3,3,6];print"C242_CASE26_BASE17_READY",Cputime(tt);
input:=eval Read("work/deg10/next_cases_20261030/C242/after_cheap_stabilizers_unified.m");row:=[r:r in input[2]|r[1] eq 26][1][2][1];assert row[1] eq 6;def:=[v[2]:v in input[1]|v[1] eq 6][1];L<a>:=ext<k5|P5!def>;ky:=row[2];b:=Evaluate(P5!ky[1..6],a);c:=Evaluate(P5!ky[7..12],a);E:=EllipticCurve([L|0,c,0,b*(1-b-c),0]);P:=E![0,0];Q:=E![b,b];S:=6*Q;R2:=21*Q;Z2:=P+R2;ai:=aInvariants(E);x:=S[1];y:=S[2];m:=(3*x^2+2*ai[2]*x+ai[4])/(2*y);a2p:=3*x+ai[2]-m^2;a3p:=2*y;kap:=a3p/a2p;b7:=-a2p^3/a3p^2;c7:=1-2*m*a2p/a3p;v7:=b7/c7;us:=4*(P[1]-R2[1])*(P[1]-Z2[1])*(R2[1]-Z2[1])/(kap^6*v7^4*(v7-1)^4);umin:=MinimalPolynomial(us,k5);assert Degree(umin) eq 2;Dq:=Divisor(Evaluate(umin,u5));pp,mm:=Support(Dq);qs5:=[pp[i]:i in [1..#pp]|mm[i] gt 0 and Degree(pp[i]) eq 2 and Valuation(j5,pp[i]) gt 0];assert #qs5 gt 0;print"C242_CASE26_Q2_CONTAINING_SET",#qs5,Coefficients(umin),Cputime(tt);
k<zeta25>:=GF(25);F,coords:=MakeField(k,model);assert Genus(F) eq 17 and Degree(F) eq 8 and Genus(BaseField(F)) eq 0;B:=BaseField(F);Pk:=PolynomialRing(k);OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);OF5:=MaximalOrderFinite(F5);OI5:=MaximalOrderInfinite(F5);
function ExtFunction(f)
 return F![B!(Pk!Coefficients(Numerator(v)))/B!(Pk!Coefficients(Denominator(v))):v in Eltseq(F5!f)];
end function;
function ExtendPlace(p)
 a,b:=TwoGenerators(p);O5:=IsFinite(p) select OF5 else OI5;assert ideal<O5|a,b> eq Ideal(p);O:=IsFinite(p) select OF else OI;I:=ideal<O|ExtFunction(a),ExtFunction(b)>;fa:=Factorization(I);assert &and[t[2] eq 1:t in fa];D:=&+[Divisor(Place(t[1])):t in fa];assert Degree(D) eq Degree(p);return D;
end function;
u:=ExtFunction(u5);s:=u^2+8;D:=s^2-3*s+9;A:=2*s^3-9*s^2+27*s-27;assert A^2-4*D^3 eq -27*D^2;Rk<Z>:=PolynomialRing(k);roots:=Roots(Z^2+Z+1);assert #roots eq 2;zeta:=roots[1][1];assert zeta^3 eq 1 and zeta ne 1;fk:=(A+3*(zeta-zeta^2)*D)/2;df:=Divisor(fk);sp,ms:=Support(df);branch:=&+[Degree(sp[i]):i in [1..#sp]|ms[i] mod 3 ne 0];assert branch eq 12;chars:=[DivisorGroup(F)!0] cat [&+[Floor(j*ms[i]/3)*Divisor(sp[i]):i in [1..#sp]]:j in [1,2]];assert [Degree(D):D in chars] eq [0,-6,-6];print"C242_CASE26_KUMMER_CHARACTERS_READY",branch,Cputime(tt);
rat5:=[p:p in cs5|Degree(p) eq 1];cub5:=[p:p in cs5|Degree(p) eq 3];assert #rat5 eq 6 and #cub5 eq 2;rat:=[ExtendPlace(p):p in rat5];cub:=[ExtendPlace(p):p in cub5];candidate:=[* *];
for i in [1..#qs5] do
 q:=ExtendPlace(qs5[i]);ps,es:=Support(q);assert #ps eq 2 and &and[Degree(p) eq 1:p in ps] and es eq [1,1];status:=[];
 for p in ps do val:=Valuation(fk,p);if val mod 3 ne 0 then Append(~status,false);continue;end if;a,b:=TwoGenerators(p);va:=Valuation(a,p);vb:=Valuation(b,p);assert Minimum(va,vb) eq 1;unif:=va eq 1 select a else b;unit:=fk/unif^val;assert Valuation(unit,p) eq 0;vv:=Evaluate(unit,p);Rr<T>:=PolynomialRing(Parent(vv));Append(~status,IsIrreducible(T^3-vv));end for;
 assert status[1] eq status[2];if &and status then Append(~candidate,[*i,q*]);end if;
end for;assert #candidate gt 0;print"C242_CASE26_INERT_Q2_CANDIDATES",#candidate,Cputime(tt);
function EC(a)return[Integers()!x:x in Eltseq(k!a,k5)];end function;
function ER(v)return[[EC(x):x in Coefficients(Numerator(v))],[EC(x):x in Coefficients(Denominator(v))]];end function;
function EF(v)return[ER(x):x in Eltseq(F!v)];end function;
function EP(p)
 a,b:=TwoGenerators(p);O:=IsFinite(p) select OF else OI;assert ideal<O|a,b> eq Ideal(p);return[*IsFinite(p) select 1 else 0,Degree(p),[EF(a),EF(b)]*];
end function;
function ED(D)
 ps,ms:=Support(D);return[*[*EP(ps[i]),ms[i]*]:i in [1..#ps]*];
end function;
function RRProof(D)
 bb:=Basis(D:Reduction:=false,Simplification:="None");Ifin,Iinf:=Ideals(D);assert Order(Ifin) eq OF and Order(Iinf) eq OI;bf:=[F!v:v in Basis(Ifin)];bi:=[F!v:v in Basis(Iinf)];assert #bf eq 8 and #bi eq 8;assert ideal<OF|bf> eq Ifin and ideal<OI|bi> eq Iinf;return[*#bb,ED(D),[EF(v):v in bf],[EF(v):v in bi],[EF(v):v in bb]*];
end function;
packet:=[*Coefficients(DefiningPolynomial(k,k5)),[ER(v):v in Coefficients(DefiningPolynomial(F))],[EF(v):v in coords],EF(u),EF(fk),[ED(D):D in chars],[ED(D):D in rat],[ED(D):D in cub],[*[*r[1],ED(r[2])*]:r in candidate*],[*def,ky,Coefficients(umin),[Integers()!x:x in Eltseq(us,k5)]*]*];Write("work/deg10/next_cases_20261030/C242/geometry/case26_section_model.m",Sprint(packet):Overwrite:=true);
records:=[* *];fail:=[* *];
for cand in candidate do qi:=cand[1];q:=cand[2];for i in [1..2],j in [1..6] do AA:=q+cub[i]+rat[j];assert Degree(AA) eq 6;proof:=[*RRProof(AA+DD):DD in chars*];dims:=[r[1]:r in proof];Append(~records,[*qi,i,j,dims,proof*]);if &+dims ne 1 then Append(~fail,[qi,i,j]);end if;print"C242_CASE26_SECTION_RR",qi,i,j,dims,Cputime(tt);end for;end for;
Write("work/deg10/next_cases_20261030/C242/geometry/case26_section_RR.m",Sprint([*records,fail*]):Overwrite:=true);printf"C242_CASE26_SECTION_COMPLETE q2=%o patterns=%o failures=%o cpu=%o\n",#candidate,#records,#fail,Cputime(tt);quit;
