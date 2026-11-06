SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(10*10^9);Alarm(500);tt:=Cputime();AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");Attach("work/deg10/continuation_20261031b/C56/inputs/fast_hecke.m");load "work/deg10/continuation_20261031b/Y228_shared/fast_primitive_inverse.m";root:="work/deg10/continuation_20261031b/C56/";src:=eval Read(root cat "final_original_source_subset.m");tag:="operator_only_s"cat IntegerToString(SHARD);mdraw:=Read(root cat "cover_v4_A5_s0_model.m");md:=eval mdraw;k:=GF(3);PK:=PolynomialRing(k);Y:=MDX11(2,28,k);CC:=Curve(Y);FF:=FunctionField(CC);F0,mp:=AlgorithmicFunctionField(FF);F:=RationalExtensionRepresentation(F0);assert Genus(F) eq 25;bc:=[F!mp(v):v in Y`_coordinates];U:=F!mp(FF.1);V:=F!mp(FF.2);OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
function ER(v)return [[Integers()!a:a in Coefficients(Numerator(v))],[Integers()!a:a in Coefficients(Denominator(v))]];end function;
function EF(v)return [ER(a):a in Eltseq(F!v)];end function;
assert [ER(v):v in Coefficients(DefiningPolynomial(F))] eq md[1] and [EF(v):v in bc] eq md[2];Write(root cat tag cat "_model.m",mdraw:Overwrite:=true);
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
   Append(~pointrows,[*id,d,ky,ps,1,pd,Q,he,Physical(E,P),Physical(E,Q*P)*]);ix:=#pointrows;cache[id]:=ix;Write(root cat tag cat "_points_progress.m",Sprint(pointrows):Overwrite:=true);
  end if;Append(~inds,ix);
 end for;
 row:=[*pointrows,cas,inds*];Write(root cat tag cat "_case"cat IntegerToString(cas[1])cat "_operator.m",Sprint(row):Overwrite:=true);Append(~out,[*cas[1],inds*]);print "C56_BASE_ONLY_OPERATOR_COMPLETE",cas[1],Cputime(tt);
end for;Write(root cat tag cat "_summary.m",Sprint(out):Overwrite:=true);print "D31B_C56_BASE_OPERATOR_BUILDER_COMPLETE",SHARD,#out,Cputime(tt);quit;
