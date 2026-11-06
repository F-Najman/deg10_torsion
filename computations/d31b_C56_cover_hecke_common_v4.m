SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(12*10^9);Alarm(1150);tt:=Cputime();AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");Attach("work/deg10/continuation_20261031b/C56/inputs/fast_hecke.m");load "work/deg10/continuation_20261031b/Y228_shared/fast_primitive_inverse.m";root:="work/deg10/continuation_20261031b/C56/";src:=eval Read(root cat "final_original_source_subset.m");tag:="cover_v4_A"cat IntegerToString(Q)cat "_s"cat IntegerToString(SHARD);
k:=GF(3);PK:=PolynomialRing(k);Y:=MDX11(2,28,k);CC:=Curve(Y);FF:=FunctionField(CC);F0,mp:=AlgorithmicFunctionField(FF);F:=RationalExtensionRepresentation(F0);assert Genus(F) eq 25;bc:=[F!mp(v):v in Y`_coordinates];U:=F!mp(FF.1);V:=F!mp(FF.2);B,C:=Explode(bc);ram:=&+[Degree(p):p in Zeros(B) cat Poles(B)|IsOdd(Valuation(B,p))];assert ram eq 24;R<Z>:=PolynomialRing(F);assert IsIrreducible(Z^2-B);L<z>:=FunctionField(Z^2-B);assert Genus(L) eq 61;EL:=EllipticCurve([L|0,-2*C,0,C^2-4*B*(1-B-C),0]);xx:=2*B+C+2*z;PP:=EL![xx,2*z*xx];assert 56*PP eq EL!0 and 28*PP ne EL!0 and 8*PP ne EL!0;tb:=TateNormalForm_bc(EL,PP);OF:=MaximalOrderFinite(L);OI:=MaximalOrderInfinite(L);print "C56_COVER_MODEL_READY",Q,SHARD,Cputime(tt);
function ER(v)return [[Integers()!a:a in Coefficients(Numerator(v))],[Integers()!a:a in Coefficients(Denominator(v))]];end function;
function EF(v)return [ER(a):a in Eltseq(F!v)];end function;
function LF(v)return [EF(a):a in Eltseq(L!v)];end function;
function KV(vals)return Y228PrimitiveJointKey(vals);end function;
function CV(v)return [Integers()!a:a in Eltseq(v,k)];end function;
function Physical(E,P)return[*[Integers()!a:a in Coefficients(DefiningPolynomial(BaseRing(E),k))],[CV(a):a in aInvariants(E)],[CV(P[1]),CV(P[2])],KV(TateNormalForm_bc(E,P))*];end function;
model:=[*[ER(a):a in Coefficients(DefiningPolynomial(F))],[EF(a):a in bc],[EF(U),EF(V)],[[EF(a):a in Coefficients(DefiningPolynomial(L))]],[LF(a):a in tb],[LF(a):a in aInvariants(EL)],[LF(PP[1]),LF(PP[2])],ram*];Write(root cat tag cat "_model.m",Sprint(model):Overwrite:=true);
function Get(E,P)
 st:=Cputime();K:=BaseRing(E);RK<T>:=PolynomialRing(K);ET,phi:=IsogenyFromKernel(E,T-(28*P)[1]);QQ:=phi(P);assert Order(QQ) eq 28;dual:=DualIsogeny(phi);ai:=aInvariants(ET);ss:=[ET![r[1],-(ai[1]*r[1]+ai[3])/2]:r in Roots(DivisionPolynomial(ET,2))];good:=[S:S in ss|dual(S) eq E!0];assert #good eq 1 and good[1] ne 14*QQ;vals:=DerickxNormalForm_bc(ET,good[1],QQ);p:=Y228FastPrimitivePlace(DefiningPolynomial(CC),F,U,V,bc,vals);gg,hh:=TwoGenerators(p);II:=ideal<OF|L!gg,L!hh>;fa:=Factorization(II);qs:=[Place(v[1]):v in fa|v[2] gt 0];assert &+[Degree(q):q in qs] eq 2*Degree(p) and &and[v[2] eq 1:v in fa|v[2] gt 0];key:=KV(TateNormalForm_bc(E,P));matched:=[q:q in qs|KV([Evaluate(f,q):f in tb]) eq key];assert #matched eq 1;out:=matched[1];assert Degree(out) eq key[1];print "C56_COVER_GET",Q,SHARD,Degree(p),Degree(out),Cputime(st),Cputime(tt);assert ideal<MaximalOrderFinite(F)|gg,hh> eq Ideal(p);desc:=[*Degree(out),Degree(p),IsFinite(p) select 1 else 0,[EF(gg),EF(hh)],key*];return out,desc;
end function;
fields:=AssociativeArray();for row in src[1] do fields[row[1]]:=ext<k|PK!row[2]>;end for;lookup:=AssociativeArray();for row in src[2] do lookup[row[1]]:=row;end for;cache:=AssociativeArray();divcache:=AssociativeArray();pointrows:=[* *];out:=[* *];
selected:=[c:c in src[3]|(Q eq 11 select c[1] eq 65 else c[1] ne 65)];
for ci in[1..#selected] do
 if (ci-1) mod NSHARDS ne SHARD then continue;end if;cas:=selected[ci];D:=DivisorGroup(L)!0;inds:=[];
 for id in cas[2] do
  if IsDefined(cache,id) then ix:=cache[id];A:=divcache[id];else
   row:=lookup[id];d:=row[2];ky:=row[3];K:=fields[d];a:=K.1;b0:=Evaluate(PK!ky[1..d],a);c0:=Evaluate(PK!ky[d+1..2*d],a);E:=EllipticCurve([K|1-c0,-b0,-b0,0,0]);P:=E![0,0];assert Order(P) eq 56;ps,psc:=Get(E,P);pd,pdc:=Get(E,Q*P);assert Degree(ps) eq d and Degree(pd) eq d;A:=-Divisor(ps)-Q*Divisor(pd);he:=[* *];
   fast:=MDIsogeniesFast(E,P,Q);assert &+[v[3]:v in fast] eq Q+1;
   for rr in fast do ET,PT,e:=Explode(rr);assert Order(PT) eq 56;p,pc:=Get(ET,PT);n:=d*e div Degree(p);assert n*Degree(p) eq d*e;A+:=n*Divisor(p);Append(~he,[*e,n,pc,Physical(ET,PT)*]);end for;
   assert Degree(A) eq 0;Append(~pointrows,[*id,d,ky,psc,1,pdc,Q,he,Physical(E,P),Physical(E,Q*P)*]);ix:=#pointrows;cache[id]:=ix;divcache[id]:=A;Write(root cat tag cat "_points_progress.m",Sprint(pointrows):Overwrite:=true);
  end if;D+:=A;Append(~inds,ix);
 end for;
 assert Degree(D) eq 0;print "C56_COVER_RR_START",Q,cas[1],Cputime(tt);Write(root cat tag cat "_operator_progress.m",Sprint([*pointrows,cas,inds*]):Overwrite:=true);st:=Cputime();h0:=#Basis(D:Reduction:=false,Simplification:="None");Append(~out,[*cas[1],cas[2],cas[3],inds,h0*]);Write(root cat tag cat "_progress.m",Sprint([*pointrows,out*]):Overwrite:=true);print "C56_COVER_CASE",Q,cas[1],h0,Cputime(st),Cputime(tt);
end for;
Write(root cat tag cat ".m",Sprint([*pointrows,out*]):Overwrite:=true);print "D31B_C56_COVER_TOWER_HECKE_COMPLETE",Q,SHARD,#out,[[v[1],v[5]]:v in out],Cputime(tt);quit;
