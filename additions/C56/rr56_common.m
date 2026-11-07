// Riemann-Roch test for the two divisors D' with label 27 of the table of Lemma "lem:sv-56-rr" (cases 492 and 737
// of continuation_20261031b/C56/final_original_source_subset.m): l(A_5 D') = 0 on X_1(56) over F_3.
// Part 1 computes A_5 D' = T_5 D' - D' - 5<5>D' from the isogenies of degree 5 on the moduli model, as in
// computations/d31b_C56_base_operator_common.m. Part 2 computes l(A_5 D') on the double cover
// X_1(56) -> X_1(2,28), k(X_1(56)) = k(X_1(2,28))(z), z^2 = b, as in computations/d31b_C56_kummer_RR_direct_v4.m
// (used for the other 9 divisors of the table). Both parts are copied from these scripts, with paths relative to the
// repository root and without writing files.
// Run from the repository root with CASE set, e.g. additions/C56/rr56_case492.m.
// ---- Part 1: the divisor A_5 D'
SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(12*10^9);Alarm(900);tt:=Cputime();
assert assigned CASE and CASE in [492,737];
AttachSpec("shared/mdmagma/v2/mdmagma.spec");Attach("continuation_20261031b/C56/inputs/fast_hecke.m");load "continuation_20261031b/Y228_shared/fast_primitive_inverse.m";root:="continuation_20261031b/C56/";src:=eval Read(root cat "final_original_source_subset.m");tag:="rr56_case"cat IntegerToString(CASE);Q:=5;CASES:=[CASE];mdraw:=Read(root cat "cover_v4_A5_s0_model.m");md:=eval mdraw;k:=GF(3);PK:=PolynomialRing(k);Y:=MDX11(2,28,k);CC:=Curve(Y);FF:=FunctionField(CC);F0,mp:=AlgorithmicFunctionField(FF);F:=RationalExtensionRepresentation(F0);assert Genus(F) eq 25;bc:=[F!mp(v):v in Y`_coordinates];U:=F!mp(FF.1);V:=F!mp(FF.2);OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
function ER(v)return [[Integers()!a:a in Coefficients(Numerator(v))],[Integers()!a:a in Coefficients(Denominator(v))]];end function;
function EF(v)return [ER(a):a in Eltseq(F!v)];end function;
assert [ER(v):v in Coefficients(DefiningPolynomial(F))] eq md[1] and [EF(v):v in bc] eq md[2];
function KV(v)return Y228PrimitiveJointKey(v);end function;
function CV(v)return [Integers()!a:a in Eltseq(v,k)];end function;
function Physical(E,P)return[*[Integers()!a:a in Coefficients(DefiningPolynomial(BaseRing(E),k))],[CV(a):a in aInvariants(E)],[CV(P[1]),CV(P[2])],KV(TateNormalForm_bc(E,P))*];end function;
function GetDesc(E,P)
 st:=Cputime();KM:=BaseRing(E);R<T>:=PolynomialRing(KM);ET,ph:=IsogenyFromKernel(E,T-(28*P)[1]);QQ:=ph(P);assert 28*QQ eq ET!0 and 14*QQ ne ET!0 and 4*QQ ne ET!0;dual:=DualIsogeny(ph);ai:=aInvariants(ET);ss:=[ET![r[1],-(ai[1]*r[1]+ai[3])/2]:r in Roots(DivisionPolynomial(ET,2))];good:=[S:S in ss|dual(S) eq E!0];assert #good eq 1 and good[1] ne 14*QQ;vals:=DerickxNormalForm_bc(ET,good[1],QQ);p:=Y228FastPrimitivePlace(DefiningPolynomial(CC),F,U,V,bc,vals);gg,hh:=TwoGenerators(p);assert ideal<OF|gg,hh> eq Ideal(p);bb,cc:=Explode(vals);key:=KV(TateNormalForm_bc(E,P));e:=Degree(p);assert key[1] in[e,2*e];EE:=EllipticCurve([KM|0,-2*cc,0,cc^2-4*bb*(1-bb-cc),0]);want:=TateNormalForm_bc(E,P);found:=[];
 for r in Roots(R![-bb,0,1]) do zr:=r[1];xx:=2*bb+cc+2*zr;PP:=EE![xx,2*zr*xx];if TateNormalForm_bc(EE,PP) eq want then Append(~found,zr);end if;end for;assert #found eq 1;z0:=found[1];assert LCM(e,Degree(MinimalPolynomial(z0,k))) eq key[1];print "C56_BASE_ONLY_PRIME",e,key[1],Cputime(st),Cputime(tt);return[*key[1],e,IsFinite(p) select 1 else 0,[EF(gg),EF(hh)],key,[CV(z0),CV(bb),CV(cc)]*];
end function;
fields:=AssociativeArray();for row in src[1] do fields[row[1]]:=ext<k|PK!row[2]>;end for;lookup:=AssociativeArray();for row in src[2] do lookup[row[1]]:=row;end for;cache:=AssociativeArray();pointrows:=[* *];out:=[* *];
for cas in src[3] do
 if not cas[1] in CASES then continue;end if;inds:=[];
 for id in cas[2] do
  if IsDefined(cache,id) then ix:=cache[id];else
   row:=lookup[id];d:=row[2];ky:=row[3];KM:=fields[d];a:=KM.1;b:=Evaluate(PK!ky[1..d],a);c:=Evaluate(PK!ky[d+1..2*d],a);E:=EllipticCurve([KM|1-c,-b,-b,0,0]);P:=E![0,0];assert 56*P eq E!0 and 28*P ne E!0 and 8*P ne E!0;ps:=GetDesc(E,P);pd:=GetDesc(E,Q*P);assert ps[1] eq d and pd[1] eq d;he:=[* *];fast:=MDIsogeniesFast(E,P,Q);assert &+[h[3]:h in fast] eq Q+1;
   for row in fast do ET,PT,e:=Explode(row);assert 56*PT eq ET!0 and 28*PT ne ET!0 and 8*PT ne ET!0;pp:=GetDesc(ET,PT);n:=d*e div pp[1];assert n*pp[1] eq d*e;Append(~he,[*e,n,pp,Physical(ET,PT)*]);end for;
   Append(~pointrows,[*id,d,ky,ps,1,pd,Q,he,Physical(E,P),Physical(E,Q*P)*]);ix:=#pointrows;cache[id]:=ix;
  end if;Append(~inds,ix);
 end for;
 row:=[*pointrows,cas,inds*];opraw:=Sprint(row);Append(~out,[*cas[1],inds*]);print "C56_BASE_ONLY_OPERATOR_COMPLETE",cas[1],Cputime(tt);
end for;assert #out eq 1 and out[1][1] eq CASE;print "RR56_OPERATOR_COMPLETE",CASE,Cputime(tt);

// ---- Part 2: the Riemann-Roch space of A_5 D'
md:=eval mdraw;op:=eval opraw;pointrows,cas,inds:=Explode(op);tag:="direct_case"cat IntegerToString(cas[1])cat "_A"cat IntegerToString(pointrows[inds[1]][7]);
k:=GF(3);PK:=PolynomialRing(k);K<t>:=FunctionField(k);R<X>:=PolynomialRing(K);
function DR(v)return K!(PK!v[1])/(K!(PK!v[2]));end function;
F<u>:=FunctionField(R![DR(v):v in md[1]]);assert Genus(F) eq 25;OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
function DF(v)return F![DR(a):a in v];end function;
function ER(v)return [[Integers()!a:a in Coefficients(Numerator(v))],[Integers()!a:a in Coefficients(Denominator(v))]];end function;
function EF(v)return [ER(a):a in Eltseq(F!v)];end function;
function PD(p)g,h:=TwoGenerators(p);OO:=IsFinite(p) select OF else OI;assert ideal<OO|g,h> eq Ideal(p);return[*Degree(p),IsFinite(p) select 1 else 0,[EF(g),EF(h)]*];end function;
function DP(v)OO:=v[3] eq 1 select OF else OI;gs:=[DF(a):a in v[4]];II:=ideal<OO|gs>;assert #gs eq 2 and IsPrime(II);p:=Place(II);assert Degree(p) eq v[2];return p;end function;
function KV(v)return Y228PrimitiveJointKey(v);end function;
function CV(v)return [Integers()!a:a in Eltseq(v,k)];end function;
function COORD(v,MQ,piv,MI)vv:=Vector(k,Eltseq(v,k));sol:=Vector(k,[vv[j]:j in piv])*MI;assert sol*MQ eq vv;return Eltseq(sol);end function;
function MM(M)return[*Nrows(M),Ncols(M),[Integers()!a:a in Eltseq(M)]*];end function;
B,C:=Explode([DF(v):v in md[2]]);assert DF(md[3][1]) eq F.1 and DF(md[3][2]) eq F!K.1;fcoeff:=Coefficients(DefiningPolynomial(F));
function VR(v,x)return Evaluate(PK!Numerator(v),x)/Evaluate(PK!Denominator(v),x);end function;
function EVALF(f,x,y)co:=Eltseq(F!f);return &+[VR(co[i],y)*x^(i-1):i in[1..#co]];end function;
function EVALPOL(x,y)return &+[VR(fcoeff[i],y)*x^(i-1):i in[1..#fcoeff]];end function;
function EVALDU(x,y)return &+[Parent(x)!(i-1)*VR(fcoeff[i],y)*x^(i-2):i in[2..#fcoeff]];end function;
function EVALDV(x,y)return &+[VR(Derivative(fcoeff[i]),y)*x^(i-1):i in[1..#fcoeff]];end function;
function PrimitiveUV(BB,CC)
 KK:=Parent(BB);rr<W>:=PolynomialRing(KK);tt:=2*BB/(1-CC)-1;qq:=((1-CC)*tt^2/2-1)/tt;if qq ne 1 then uu:=W;vv:=((qq+1)*W+2)/(qq-1);else uu:=rr!(-1);vv:=W;end if;nt:=uu^3+uu^2*vv+2*uu^2-uu*vv^2-vv^3-2*vv^2;dt:=uu^3+uu^2*vv+2*uu^2+uu*vv^2+2*uu*vv+vv^3+2*vv^2-2;sol:=[];
 for root in Roots(nt-tt*dt) do aa:=Evaluate(uu,root[1]);bb:=Evaluate(vv,root[1]);if aa eq bb then continue;end if;den:=aa^3+aa^2*bb+2*aa^2+aa*bb^2+2*aa*bb+bb^3+2*bb^2-2;if den eq 0 then continue;end if;t0:=(aa^3+aa^2*bb+2*aa^2-aa*bb^2-bb^3-2*bb^2)/den;q0:=(-aa-bb-2)/(aa-bb);if EVALPOL(aa,bb) eq 0 and t0 eq tt and q0 eq qq then Append(~sol,[aa,bb]);end if;end for;assert #sol eq 1;return sol[1];
end function;
places:=[];groups:=[* *];
procedure AddTerm(~places,~groups,desc,n,physical)
 p:=DP(desc);ii:=Index(places,p);if ii eq 0 then Append(~places,p);Append(~groups,[*[],[],[* *],[* *]*]);ii:=#places;end if;gr:=groups[ii];key:=Sprint(desc[5]);jj:=Index(gr[1],key);
 if jj eq 0 then Append(~gr[1],key);Append(~gr[2],n);Append(~gr[3],desc);Append(~gr[4],physical);else gr[2][jj]+:=n;end if;groups[ii]:=gr;
end procedure;
for ix in inds do row:=pointrows[ix];Q:=row[7];AddTerm(~places,~groups,row[4],-1,row[9]);AddTerm(~places,~groups,row[6],-Q,row[10]);for h in row[8] do AddTerm(~places,~groups,h[3],h[2],h[4]);end for;end for;
A0:=DivisorGroup(F)!0;ac:=[];ddeg:=0;
for i in[1..#places] do a:=Max([0]cat groups[i][2]);Append(~ac,a);A0+:=a*Divisor(places[i]);for j in[1..#groups[i][2]] do ddeg+:=groups[i][2][j]*groups[i][3][j][1];end for;end for;assert ddeg eq 0;
psB:=Zeros(B) cat Poles(B);Rhalf:=&+[Integers()!Floor(Rationals()!Valuation(B,p)/2)*Divisor(p):p in psB];assert Degree(Rhalf) eq -12;V0,v0:=RiemannRochSpace(A0);V1,v1:=RiemannRochSpace(A0+Rhalf);f0:=[v0(w):w in Basis(V0)];f1:=[v1(w):w in Basis(V1)];n0:=#f0;n1:=#f1;nb:=n0+n1;print "C56_KUMMER_BASE_SPACES",cas[1],Degree(A0),n0,n1,Cputime(tt);
ALL:=ZeroMatrix(k,nb,0);tables:=[* *];
for gi in[1..#places] do
 p:=places[gi];e:=Degree(p);a:=ac[gi];gr:=groups[gi];if a eq 0 and &and[n eq 0:n in gr[2]] then continue;end if;
 physical:=gr[4][1];PF:=PolynomialRing(k);KM<aa>:=ext<k|PF!physical[1]>;decode:=func<v|Evaluate(PF!v,aa)>;E:=EllipticCurve([decode(v):v in physical[2]]);P:=E![decode(v):v in physical[3]];assert 56*P eq E!0 and 28*P ne E!0 and 8*P ne E!0;RK<T>:=PolynomialRing(KM);ET,ph:=IsogenyFromKernel(E,T-(28*P)[1]);QQ:=ph(P);dual:=DualIsogeny(ph);ai:=aInvariants(ET);ss:=[ET![r[1],-(ai[1]*r[1]+ai[3])/2]:r in Roots(DivisionPolynomial(ET,2))];S:=[s:s in ss|dual(s) eq E!0][1];bcphys:=DerickxNormalForm_bc(ET,S,QQ);bb,cc:=Explode(bcphys);EE:=EllipticCurve([KM|0,-2*cc,0,cc^2-4*bb*(1-bb-cc),0]);want:=TateNormalForm_bc(E,P);roots:=Roots(RK![-bb,0,1]);assert #roots eq 2;z0s:=[];
 for rr in roots do zr:=rr[1];xx:=2*bb+cc+2*zr;PP:=EE![xx,2*zr*xx];if TateNormalForm_bc(EE,PP) eq want then Append(~z0s,zr);end if;end for;assert #z0s eq 1;z0:=z0s[1];dq:=gr[3][1][1];assert dq in[e,2*e];
 keys:=[* *];coeffs:=[];qdegrees:=[];qsigns:=[];
 for sg in (dq eq e select [1,-1] else [1]) do zr:=sg*z0;xx:=2*bb+cc+2*zr;PP:=EE![xx,2*zr*xx];key:=KV(TateNormalForm_bc(EE,PP));assert key[1] eq dq;idx:=Index(gr[1],Sprint(key));n:=idx eq 0 select 0 else gr[2][idx];Append(~keys,key);Append(~coeffs,n);Append(~qdegrees,dq);Append(~qsigns,sg);end for;assert &and[key in [Sprint(z):z in keys]:key in gr[1]];needs:=[a-n:n in coeffs];assert &and[n ge 0:n in needs];maxneed:=Max(needs);if maxneed eq 0 then continue;end if;
 uv0:=PrimitiveUV(bb,cc);uu0,vv0:=Explode(uv0);du0:=EVALDU(uu0,vv0);dv0:=EVALDV(uu0,vv0);assert du0 ne 0 or dv0 ne 0;prec:=128;LS<uni>:=LaurentSeriesRing(KM,prec);UU:=LS!uu0+O(uni^prec);VV:=LS!vv0+O(uni^prec);which:=du0 ne 0 select 2 else 1;if which eq 2 then VV+:=uni;else UU+:=uni;end if;
 for it in[1..8] do ff:=EVALPOL(UU,VV);if which eq 2 then dd:=EVALDU(UU,VV);assert Valuation(dd) eq 0;UU-:=ff/dd;else dd:=EVALDV(UU,VV);assert Valuation(dd) eq 0;VV-:=ff/dd;end if;end for;assert Valuation(EVALPOL(UU,VV)) ge prec-2;assert Coefficient(UU,0) eq uu0 and Coefficient(VV,0) eq vv0;
 assert &and[Valuation(EVALF(DF(v),UU,VV)) gt 0:v in gr[3][1][4]];sB:=EVALF(B,UU,VV);assert Valuation(sB) eq 0 and Coefficient(sB,0) eq bb;BBs:=[Coefficient(sB,n):n in[0..maxneed-1]];zz:=[z0];
 for n in[1..maxneed-1] do cross:=KM!0;for i in[1..n-1] do cross+:=zz[i+1]*zz[n-i+1];end for;Append(~zz,(BBs[n+1]-cross)/(2*z0));end for;assert &and[&+[zz[i+1]*zz[n-i+1]:i in[0..n]] eq BBs[n+1]:n in[0..maxneed-1]];
 s0:=[uni^a*EVALF(f,UU,VV):f in f0];s1:=[uni^a*EVALF(f,UU,VV):f in f1];assert &and[Valuation(ss) ge 0 and AbsolutePrecision(ss) ge maxneed:ss in s0 cat s1];H0:=[[Coefficient(ss,n):n in[0..maxneed-1]]:ss in s0];H1:=[[Coefficient(ss,n):n in[0..maxneed-1]]:ss in s1];ZH1:=[[&+[zz[i+1]*row[n-i+1]:i in[0..n]]:n in[0..maxneed-1]]:row in H1];blocks:=[* *];
 for qi in[1..#keys] do need:=needs[qi];if need eq 0 then continue;end if;sg:=qsigns[qi];zr:=sg*z0;xx:=2*bb+cc+2*zr;PP:=EE![xx,2*zr*xx];uv:=TateNormalForm_bc(EE,PP);bits:=keys[qi][2];gq:=&+[KM!bits[i]*uv[i]:i in[1,2]];polq:=MinimalPolynomial(gq,k);assert Degree(polq) eq dq;LQ<vq>:=ext<k|polq>;hq:=hom<LQ->KM|gq>;MQ:=Matrix(k,[Eltseq(gq^i,k):i in[0..dq-1]]);assert Rank(MQ) eq dq;EQ:=EchelonForm(MQ);piv:=[Min([j:j in[1..Ncols(EQ)]|EQ[i,j] ne 0]):i in[1..dq]];MI:=Matrix(k,[[MQ[i,j]:j in piv]:i in[1..dq]])^-1;matrows:=[];
  for row in H0 do vals:=[];for n in[1..need] do vals cat:=COORD(row[n],MQ,piv,MI);end for;Append(~matrows,vals);end for;
  for row in ZH1 do vals:=[];for n in[1..need] do value:=sg*row[n];vals cat:=COORD(value,MQ,piv,MI);end for;Append(~matrows,vals);end for;
  MM0:=Matrix(k,matrows);assert Nrows(MM0) eq nb and Ncols(MM0) eq dq*need;ALL:=HorizontalJoin(ALL,MM0);Append(~blocks,[*keys[qi],coeffs[qi],sg,dq,need,CV(gq),[Integers()!x:x in Coefficients(polq)],MM(MM0)*]);
 end for;
 Append(~tables,[*PD(p),a,e,[Integers()!x:x in Coefficients(DefiningPolynomial(KM,k))],[CV(v):v in bcphys],CV(z0),[CV(v):v in BBs],[[CV(v):v in row]:row in H0],[[CV(v):v in row]:row in H1],blocks,[*which,prec,[CV(Coefficient(UU,j)):j in[0..prec-1]],[CV(Coefficient(VV,j)):j in[0..prec-1]]*]*]);print "C56_KUMMER_LOCAL_BLOCK",cas[1],gi,e,a,needs,Ncols(ALL),Cputime(tt);
end for;
assert Ncols(ALL) eq 2*Degree(A0);rk:=Rank(ALL);h0:=nb-rk;print "D31B_C56_DIRECT_KUMMER_RR_COMPLETE",cas[1],pointrows[inds[1]][7],Degree(A0),n0,n1,Ncols(ALL),rk,h0,Cputime(tt);assert cas[1] eq CASE and pointrows[inds[1]][7] eq 5 and h0 eq 0;
print "RR56_LABEL27_DONE",CASE,"deg(H) =",Degree(A0),"l(H) =",n0,"l(H+R_b) =",n1,"conditions =",Ncols(ALL),"rank =",rk,"l(A_5 D) =",h0,Cputime(tt);
quit;
