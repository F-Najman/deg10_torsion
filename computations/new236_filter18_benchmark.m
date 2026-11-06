TAG:="benchmark";SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(5*10^9);tt:=Cputime();k:=GF(5);PP:=PolynomialRing(k);
a:=eval Read("work/deg10/C236_resolution_20261030/enumeration/target18_proven_picard.m");B<t>:=FunctionField(k);RF<Y>:=PolynomialRing(B);
load "work/deg10/C236_resolution_20261030/enumeration/deps/pair_af_place_linear.m";
load "work/deg10/C236_resolution_20261030/enumeration/deps/portable_forward.m";
function Rat(v) return B!(PP!v[1])/B!(PP!v[2]);end function;
f:=RF![Rat(v):v in a[1]];f/:=LeadingCoefficient(f);F<y>:=FunctionField(f);assert Genus(F) eq 7;
function GetFunction(data) return F![Rat(v):v in data];end function;
OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);fb:=[];
for j in [1..#a[4]] do OO:=a[2][j] eq 1 select OF else OI;I:=ideal<OO|[GetFunction(v):v in a[4][j]]>;pl:=Place(I);assert Degree(pl) eq a[3][j];Append(~fb,pl);end for;
Z:=FreeAbelianGroup(#a[6][1]);G,q:=quo<Z|[Z!v:v in a[6]]>;fc:=[q(Z!v):v in a[5]];A:=Divisor(fb[a[8]]);coordsBC:=[GetFunction(v):v in a[7]];b,c:=Explode(coordsBC);tf:=2*b/(1-c)-1;qf:=((1-c)*tf^2-2)/(2*tf);uf:=1/qf;vf:=tf*(uf^2-1)/(tf*(uf^2+3)+4*uf);coordsUV:=[uf,vf];
phi,witness:=C240PortableDirectRR(F,fb,fc,G,A,Sprintf("work/deg10/C236_resolution_20261030/enumeration/filter18_%o_RR.log",TAG));hg:=[q(Z!v):v in a[9]];H:=sub<G|hg>;assert Invariants(H) eq [2,42,126];QH,qh:=quo<G|H>;assert #TorsionSubgroup(QH) eq 12;
function Key(vals,D)
 keys:=[];for i in [1..D] do Append(~keys,&cat[[Integers()!a:a in Eltseq(z,k)]:z in vals]);vals:=[z^5:z in vals];end for;return Sort(keys)[1];
end function;
function BC(E,P,Q)
 ai:=aInvariants(E);xx:=Q[1]-P[1];yy:=Q[2];return [xx^3/yy^2,(ai[2]+3*P[1])*xx^2/yy^2];
end function;
function Pair(vals)
 b,c:=Explode(vals);if c eq 1 then return false,vals;end if;tt:=2*b/(1-c)-1;if tt eq 0 then return false,vals;end if;qq:=((1-c)*tt^2-2)/(2*tt);if qq eq 0 then return false,vals;end if;uu:=1/qq;den:=tt*(uu^2+3)+4*uu;if den eq 0 then return false,vals;end if;vv:=tt*(uu^2-1)/den;return true,[uu,vv];
end function;
function ER(v) return [[Integers()!c:c in Coefficients(Numerator(v))],[Integers()!c:c in Coefficients(Denominator(v))]];end function;
function EF(v) return [ER(c):c in Eltseq(F!v)];end function;
Gn,fromn,ton:=ClassGroup(F:ReductionDivisor:=A,Proof:=true);print "C236_NATIVE_MAP_SETUP",Cputime(tt);cn:=[ton(Divisor(p)):p in fb];V:=FreeAbelianGroup(#fb);oldm:=hom<V->G|fc>;newm:=hom<V->Gn|cn>;assert Kernel(oldm) eq Kernel(newm);isom:=hom<G->Gn|[newm(G.i@@oldm):i in [1..Ngens(G)]]>;assert #Kernel(isom) eq 1;
src:=eval Read("work/deg10/C236_resolution_20261030/enumeration/source6.m");K<zz>:=ext<k|PP!src[2]>;native:=0.;direct:=0.;tested:=0;
for rec in src[9] do
 if rec[1] ne 6 then continue;end if;ky:=rec[2];b:=Evaluate(PP!ky[1..6],zz);c:=Evaluate(PP!ky[7..12],zz);E:=EllipticCurve([K|0,c,0,b*(1-b-c),0]);P:=E![0,0];Q:=E![b,b];bc:=BC(E,P,2*Q);fast,uv:=Pair(bc);pl:=C48PairAFPlaceLinear(fast select coordsUV else coordsBC,uv);e:=Degree(pl);D:=Divisor(pl)-e*A;t:=Cputime();v:=ton(D);native+:=Cputime(t);t:=Cputime();u:=phi(D);direct+:=Cputime(t);assert v eq isom(u);tested+:=1;if tested eq 30 then break;end if;
end for;fprintf witness,"]\n";Flush(witness);print "C236_NATIVE_MAP_BENCHMARK_COMPLETE",tested,native,direct,Cputime(tt);quit;
