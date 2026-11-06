SetColumns(0);SetNthreads(1);SetQuitOnError(true);SetMemoryLimit(10*10^9);Alarm(1600);tt:=Cputime();AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");Attach("work/deg10/continuation_20261031b/Y228_shared/fast_hecke_11.m");load "work/deg10/continuation_20261031b/Y228_shared/pair_af_place_linear.m";
root:="work/deg10/continuation_20261031b/Y228_shared/";load "work/deg10/continuation_20261031b/Y228_shared/union_cases.m";k:=GF(3);PK:=PolynomialRing(k);X:=MDX11(2,28,k);CC:=Curve(X);F0,mp:=AlgorithmicFunctionField(FunctionField(CC));F:=RationalExtensionRepresentation(F0);assert Genus(F)eq 25;coords:=[F!mp(v):v in X`_coordinates];OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);tag:=IntegerToString(SHARD);print "Y228_SHARED_AF_READY",SHARD,Degree(F),Genus(F),Cputime(tt);
function ER(v)return [[Integers()!a:a in Coefficients(Numerator(v))],[Integers()!a:a in Coefficients(Denominator(v))]];end function;
function EF(v)return [ER(a):a in Eltseq(F!v)];end function;
function PD(p)a,b:=TwoGenerators(p);O:=IsFinite(p)select OF else OI;assert ideal<O|a,b>eq Ideal(p);return[*Degree(p),IsFinite(p)select 1 else 0,[EF(a),EF(b)]*];end function;
function Get(E,P,Q)return C48PairAFPlaceLinear(coords,DerickxNormalForm_bc(E,P,Q));end function;
cs:=Cusps(X);cub:=[c:c in cs|Degree(c)eq 3];assert #cub eq 4;cp:=[];for c in cub do dv:=DivisorGroup(F)!FunctionFieldDivisor(Divisor(c));ps,es:=Support(dv);assert #ps eq 1 and es eq [1]and Degree(ps[1])eq 3;Append(~cp,ps[1]);end for;
b,c:=Explode(coords);d:=b*(1-b-c);EE:=EllipticCurve([F|0,c,0,d,0]);QQ:=EE![b,b];fs:=[b,d,c^2-4*d]cat[(j*QQ)[1]:j in[1..14]];sig:=[[Valuation(v,p):v in fs]:p in cp];ss:=Sort(Setseq(Seqset(sig)));assert #ss eq 2;pairs:=[[i:i in[1..4]|sig[i]eq s]:s in ss];assert &and[#p eq 2:p in pairs];vs:=[];for pair in pairs do for eps in[-1,1]do v:=[0:i in[1..4]];v[pair[1]]:=-eps;v[pair[2]]:=eps;Append(~vs,v);end for;end for;corrections:=Sort(Setseq(Seqset([[0,0,0,0]]cat[[v[j]+5*w[j]:j in[1..4]]:v,w in vs])));assert #corrections eq 17;
Write(root cat "model_s"cat tag cat".m",Sprint([*[ER(v):v in Coefficients(DefiningPolynomial(F))],[EF(v):v in coords],[PD(p):p in cp],sig,pairs,corrections*]):Overwrite:=true);
fields:=AssociativeArray();for row in UnionFields do fields[row[1]]:=ext<k|PK!row[2]>;end for;cache:=AssociativeArray();allout:=[* *];
for ci in[1..#UnionCases]do
 if(ci-1)mod NSHARDS ne SHARD then continue;end if;id,cusp,parts:=Explode(UnionCases[ci]);assert id eq ci;D:=DivisorGroup(F)!0;obs:=[* *];
 for part in parts do
  e,mul,ky:=Explode(part);key:=[e]cat ky;
  if IsDefined(cache,key)then val:=cache[key];else
   K:=fields[e];a:=K.1;bb:=Evaluate(PK!ky[1..e],a);cc:=Evaluate(PK!ky[e+1..2*e],a);E:=EllipticCurve([K|0,cc,0,bb*(1-bb-cc),0]);P2:=E![0,0];Q:=E![bb,bb];assert Order(Q)eq 28 and P2 ne 14*Q;ps:=Get(E,P2,Q);pd:=Get(E,P2,5*Q);assert Degree(ps)eq e and Degree(pd)eq e;DD:=-Divisor(ps)-5*Divisor(pd);hec:=[* *];fast:=MDIsogeniesFastPoints(E,[P2,Q],5);assert &+[h[3]:h in fast]eq 6;
   for row in fast do ET,pts,ell:=Explode(row);assert 2*pts[1]eq ET!0 and pts[1]ne ET!0 and Order(pts[2])eq 28 and pts[1]ne 14*pts[2];T:=Get(ET,pts[1],pts[2]);n:=e*ell div Degree(T);assert n*Degree(T)eq e*ell;DD+:=n*Divisor(T);Append(~hec,[*ell,n,PD(T)*]);end for;
   assert Degree(DD)eq 0;val:=[*DD,PD(ps),PD(pd),hec*];cache[key]:=val;
  end if;D+:=mul*val[1];Append(~obs,[*e,mul,ky,val[2],val[3],val[4]*]);
 end for;assert Degree(D)eq 0;cslist:=cusp ge 3 select corrections else [[0,0,0,0]];hs:=[];
 for v in cslist do A:=D+&+[v[j]*Divisor(cp[j]):j in[1..4]];assert Degree(A)eq 0;Append(~hs,#Basis(A:Reduction:=false,Simplification:="None"));end for;
 Append(~allout,[*id,cusp,obs,cslist,hs*]);Write(root cat "A5_s"cat tag cat"_progress.m",Sprint(allout):Overwrite:=true);print "D31B_Y228_A5_CASE",SHARD,id,cusp,hs,Cputime(tt);
end for;
Write(root cat "A5_s"cat tag cat".m",Sprint(allout):Overwrite:=true);print "D31B_Y228_A5_COMPLETE",SHARD,#allout,[[v[1],v[2]]:v in allout|&+v[5]ne 0],Cputime(tt);quit;
