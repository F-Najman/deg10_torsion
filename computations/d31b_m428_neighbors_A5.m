SetColumns(0);SetNthreads(1);SetQuitOnError(true);SetMemoryLimit(10*10^9);Alarm(700);tt:=Cputime();AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");Attach("work/deg10/continuation_20261031b/Y228_shared/fast_hecke_11.m");load "work/deg10/continuation_20261031b/Y228_shared/fast_primitive_inverse.m";
root:="work/deg10/continuation_20261031b/Y228_shared/";load "work/deg10/continuation_20261031b/Y228_shared/union_cases.m";k:=GF(3);PK:=PolynomialRing(k);X:=MDX11(2,28,k);CC:=Curve(X);FF:=FunctionField(CC);F0,mp:=AlgorithmicFunctionField(FF);F:=RationalExtensionRepresentation(F0);assert Genus(F)eq 25;coords:=[F!mp(v):v in X`_coordinates];OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);tag:="m428_neighbors";print "Y228_SHARED_AF_READY",0,Degree(F),Genus(F),Cputime(tt);
function ER(v)return [[Integers()!a:a in Coefficients(Numerator(v))],[Integers()!a:a in Coefficients(Denominator(v))]];end function;
function EF(v)return [ER(a):a in Eltseq(F!v)];end function;
function PD(p)a,b:=TwoGenerators(p);O:=IsFinite(p)select OF else OI;assert ideal<O|a,b>eq Ideal(p);return[*Degree(p),IsFinite(p)select 1 else 0,[EF(a),EF(b)]*];end function;
function Get(E,P,Q)return Y228FastPrimitivePlace(DefiningPolynomial(CC),F,F!mp(FF.1),F!mp(FF.2),coords,DerickxNormalForm_bc(E,P,Q));end function;
bb0,cc0:=Explode(coords);E0:=EllipticCurve([F|0,cc0,0,bb0*(1-bb0-cc0),0]);jf:=jInvariant(E0);csAF:=Poles(jf);assert Degree(jf)eq 576 and {*Degree(p):p in csAF*}eq {*1^^12,2^^6,3^^4,6^^2*};cp:=[p:p in csAF|Degree(p)eq 3];assert #cp eq 4;print "Y228_AF_CUBIC_CUSPS_READY",0,Cputime(tt);
b,c:=Explode(coords);d:=b*(1-b-c);EE:=EllipticCurve([F|0,c,0,d,0]);QQ:=EE![b,b];fs:=[b,d,c^2-4*d];sig:=[[Valuation(v,p):v in fs]:p in cp];ss:=Sort(Setseq(Seqset(sig)));assert #ss eq 2;pairs:=[[i:i in[1..4]|sig[i]eq s]:s in ss];assert &and[#p eq 2:p in pairs];vs:=[];for pair in pairs do for eps in[-1,1]do v:=[0:i in[1..4]];v[pair[1]]:=-eps;v[pair[2]]:=eps;Append(~vs,v);end for;end for;corrections:=Sort(Setseq(Seqset([[0,0,0,0]]cat[[v[j]+5*w[j]:j in[1..4]]:v,w in vs])));assert #corrections eq 17;
Write(root cat "model_"cat tag cat".m",Sprint([*[ER(v):v in Coefficients(DefiningPolynomial(F))],[EF(v):v in coords],[PD(p):p in cp],sig,pairs,corrections*]):Overwrite:=true);
fields:=AssociativeArray();for row in UnionFields do fields[row[1]]:=ext<k|PK!row[2]>;end for;
allout:=[* *];
for ci in[27,42,43,45]do
 id,cusp,parts:=Explode(UnionCases[ci]);assert cusp eq 0 and #parts eq 1 and parts[1][1]eq 10 and parts[1][2]eq 1;ky:=parts[1][3];K:=fields[10];a:=K.1;bb:=Evaluate(PK!ky[1..10],a);cc:=Evaluate(PK!ky[11..20],a);E:=EllipticCurve([K|0,cc,0,bb*(1-bb-cc),0]);P2:=E![0,0];Q:=E![bb,bb];assert Order(Q)eq 28;Us:=DivisionPoints(P2,2);assert #Us eq 4;U:=Us[1];V:=7*Q;assert 2*U eq P2 and 2*V eq 14*Q;
 for ki in[1..3]do
 S:=ki eq 1 select P2 else (ki eq 2 select 14*Q else P2+14*Q);ET,phi:=IsogenyFromKernel(E,PolynomialRing(K).1-S[1]);
 if ki eq 1 then R2:=phi(U);T:=phi(Q);elif ki eq 2 then R2:=phi(V);T:=phi(U+8*Q);else R2:=phi(U+V);T:=phi(Q);end if;
 assert Order(T)eq 28 and R2 ne ET!0 and 2*R2 eq ET!0 and R2 ne 14*T;ps:=Get(ET,R2,T);pd:=Get(ET,R2,5*T);ns:=10 div Degree(ps);nd:=50 div Degree(pd);assert ns*Degree(ps)eq 10 and nd*Degree(pd)eq 50;D:=-ns*Divisor(ps)-nd*Divisor(pd);hec:=[* *];fast:=MDIsogeniesFastPoints(ET,[R2,T],5);assert &+[h[3]:h in fast]eq 6;
 for row in fast do ETT,pts,ell:=Explode(row);P:=Get(ETT,pts[1],pts[2]);nn:=10*ell div Degree(P);assert nn*Degree(P)eq 10*ell;D+:=nn*Divisor(P);Append(~hec,[*ell,nn,PD(P)*]);end for;
 assert Degree(D)eq 0;h0:=#Basis(D);Append(~allout,[*ci,ki,[[Integers()!v:v in Eltseq(z,k)]:z in Eltseq(U)],[[Integers()!v:v in Eltseq(z,k)]:z in aInvariants(ET)],[[Integers()!v:v in Eltseq(z,k)]:z in Eltseq(R2)],[[Integers()!v:v in Eltseq(z,k)]:z in Eltseq(T)],PD(ps),ns,PD(pd),nd,hec,h0*]);Write(root cat "m428_neighbor_A5_progress.m",Sprint(allout):Overwrite:=true);print "M428_NEIGHBOR_A5",ci,ki,h0,Cputime(tt);
 end for;
end for;
Write(root cat "m428_neighbor_A5.m",Sprint(allout):Overwrite:=true);print "D31B_M428_NEIGHBOR_A5_COMPLETE",#allout,[*[*r[1],r[2],r[12]*]:r in allout*],Cputime(tt);quit;
