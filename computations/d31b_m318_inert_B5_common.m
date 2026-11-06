SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(8*10^9);Alarm(1500);tt:=Cputime();root:="work/deg10/continuation_20261031b/M3N18/";AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");Attach("work/deg10/continuation_20261029/deps/fast_hecke_11.m");Attach(root cat "inputs/fast_full_aut_hash_m3_18_root.m");Attach(root cat "isogenies_witness.m");
k:=GF(13);X:=MDX11(3,18,k);C:=Curve(X);FF:=FunctionField(C);F,mp:=AlgorithmicFunctionField(FF);assert Genus(F)eq 10;B:=BaseField(F);O:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
function ER(f)return [[Integers()!c:c in Coefficients(Numerator(f))],[Integers()!c:c in Coefficients(Denominator(f))]];end function;
function EF(f)return [ER(c):c in Eltseq(F!f)];end function;
function PD(P)a,b:=TwoGenerators(P);assert ideal<Order(Ideal(P))|a,b>eq Ideal(P);return[*Degree(P),IsFinite(P)select 1 else 0,[EF(a),EF(b)]*];end function;
function PlacePair(qt)cp:=MDPlaceCommonZeros(X`_coordinates,qt);ps,es:=Support(FunctionFieldDivisor(Divisor(cp)));assert #ps eq 1 and es eq[1];P:=ps[1];K:=Parent(qt[1]);assert M3EighteenRootKey([Evaluate(f,cp):f in X`_coordinates],K,13,Degree(K))eq M3EighteenRootKey(qt,K,13,Degree(K));return P;end function;
function FX(vals)L:=Parent(vals[1]);theta:=PrimitiveElement(L);d:=Degree(L,k);M:=Matrix(k,[Eltseq(theta^i,k):i in[0..d-1]]);assert Rank(M)eq d;co:=[Eltseq(Solution(M,Vector(k,Eltseq(v,k)))):v in vals];return[*[Integers()!a:a in Coefficients(MinimalPolynomial(theta,k))],[[Integers()!a:a in c]:c in co]*];end function;
inputs:=eval Read(root cat "paired_final_Hecke_inputs_v2.m");out:=[* *];q:=5;RK:=PolynomialRing(k);ids:=[10,12,24,29,30,32,33,38,40,43];
for row in inputs do
 if row[1]notin ids or(Index(ids,row[1])-1)mod 3 ne SHARD then continue;end if;D:=DivisorGroup(F)!0;rs:=[* *];
 for inp in row[3]do
  d:=inp[1];K<z>:=ext<k|RK!inp[2]>;qt:=[Evaluate(RK!v,z):v in inp[3]];q0,t0:=Explode(qt);E:=EllipticCurve([q0*t0-q0+t0+2,K!0,q0*t0^2-q0*t0+t0,K!0,K!0]);P:=E![0,0];Q:=E![-t0,t0^2];assert 3*P eq E!0 and P ne E!0 and 18*Q eq E!0 and 9*Q ne E!0 and 6*Q ne E!0 and WeilPairing(P,6*Q,3)ne 1;N:=PlacePair(qt);assert Degree(N)eq d;src:=FX(qt cat aInvariants(E)cat[P[1],P[2],Q[1],Q[2]]);
  dq:=M3EighteenRootNormalizedPair(E,5*P,5*Q);dqq:=M3EighteenRootNormalizedPair(E,25*P,25*Q);Nd:=PlacePair(dq);Ndd:=PlacePair(dqq);assert Degree(Nd)eq d and Degree(Ndd)eq d;D-:=Divisor(N)+10*Divisor(Nd)+25*Divisor(Ndd);chains:=[* *];ng:=0;
  for first in MD318IsogeniesWitness(E,[P,Q],5)do
   E1,pt1,e1,EL1,ker1:=Explode(first);P1:=2*pt1[1];Q1:=pt1[2];assert WeilPairing(P1,6*Q1,3)eq BaseRing(E1)!WeilPairing(P,6*Q,3);vv1:=M3EighteenRootNormalizedPair(E1,P1,Q1);L1:=BaseRing(E1);wv1:=vv1 cat aInvariants(EL1)cat[L1!P[1],L1!P[2],L1!Q[1],L1!Q[2]]cat aInvariants(E1)cat[P1[1],P1[2],Q1[1],Q1[2]]cat Coefficients(ker1);secondrows:=[* *];nn:=0;
   for second in MD318IsogeniesWitness(E1,[P1,Q1],5)do
    E2,pt2,e2,EL2,ker2:=Explode(second);P2:=2*pt2[1];Q2:=pt2[2];assert WeilPairing(P2,6*Q2,3)eq BaseRing(E2)!WeilPairing(P1,6*Q1,3);vv2:=M3EighteenRootNormalizedPair(E2,P2,Q2);Nt:=PlacePair(vv2);mult:=d*e1*e2 div Degree(Nt);assert mult*Degree(Nt)eq d*e1*e2;D+:=mult*Divisor(Nt);L2:=BaseRing(E2);wv2:=vv2 cat aInvariants(EL2)cat[L2!P1[1],L2!P1[2],L2!Q1[1],L2!Q1[2]]cat aInvariants(E2)cat[P2[1],P2[2],Q2[1],Q2[2]]cat Coefficients(ker2);Append(~secondrows,[*PD(Nt),mult,e2,FX(wv2)*]);nn+:=e2;ng+:=e1*e2;
   end for;assert nn eq 6;Append(~chains,[*e1,FX(wv1),secondrows*]);
  end for;assert ng eq 36;Append(~rs,[*PD(N),src,PD(Nd),FX(dq),PD(Ndd),FX(dqq),chains*]);
 end for;assert Degree(D)eq 0;print "M318_INERT_B5_RR_START",row[1],Cputime(tt);h0:=#Basis(D:Reduction:=false,Simplification:="None");Append(~out,[*row[1],row[2],rs,h0*]);Write(root cat Sprintf("inert_B5_s%o_progress.m",SHARD),Sprint(out):Overwrite:=true);print "M318_INERT_B5_CASE",row[1],h0,Cputime(tt);
end for;
Write(root cat Sprintf("inert_B5_s%o.m",SHARD),Sprint([*[ER(c):c in Coefficients(DefiningPolynomial(F))],[EF(mp(f)):f in X`_coordinates],out*]):Overwrite:=true);print "D31B_M318_INERT_B5_COMPLETE",SHARD,#out,[r[1]:r in out|r[4]ne 0],Cputime(tt);quit;
