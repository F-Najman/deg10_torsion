SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(4*10^9);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
a:=eval Read("work/deg10/C232_resolution_20261030/hecke11_geometry_v2.m");k:=GF(3);PP:=PolynomialRing(k);B<t>:=FunctionField(k);RF<Y>:=PolynomialRing(B);
function Rat(v) return B!(PP!v[1])/B!(PP!v[2]);end function;
f:=RF![Rat(v):v in a[1]];f/:=LeadingCoefficient(f);F<y>:=FunctionField(f);assert Genus(F) eq 37;
function GFUN(v) return F![Rat(w):w in v];end function;
OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
function ER(v) return [[Integers()!c:c in Coefficients(Numerator(v))],[Integers()!c:c in Coefficients(Denominator(v))]];end function;
function EF(v) return [ER(c):c in Eltseq(F!v)];end function;
function EP(pl) I:=Ideal(pl);return [*Order(I) eq OF select 1 else 0,Degree(pl),[EF(v):v in Basis(I)]*];end function;
function DP(s) OO:=s[1] eq 1 select OF else OI;I:=ideal<OO|[GFUN(v):v in s[3]]>;pl:=Place(I);assert Degree(pl) eq s[2];return pl;end function;
geom:=eval Read("work/deg10/C232_resolution_20261030/exponent_cusp_geometry.m");ps:=[DP(s):s in geom[3]] cat [DP(s):s in geom[2]];
/* The two degree-eight cusps are irrelevant to a degree-four complement.
   Global Galois orbits among these 26 finite places number 23:
   16 rational, four degree-two, two split quartic, one split octic. */
coords:=[GFUN(v):v in a[2]];b,c:=Explode(coords);E:=EllipticCurve([F|0,c,0,b*(1-b-c),0]);P:=E![0,0];Q:=E![b,b];delta:=DerickxNormalForm_bc(E,P,5*Q);bb,cc:=Explode(delta);EG:=EllipticCurve([F|0,cc,0,bb*(1-bb-cc),0]);QG:=EG![bb,bb];
us:=[b,c,1-b-c,c^2-4*b*(1-b-c)];ugs:=[bb,cc,1-bb-cc,cc^2-4*bb*(1-bb-cc)];
for j in [1..9] do
 if j gt 1 then V:=j*Q;VG:=j*QG;us cat:=[V[1],V[2]];ugs cat:=[VG[1],VG[2]];end if;
 sigs:=[[Valuation(u,p):u in us]:p in ps];n:=#Seqset(sigs);print "C232_CUSP_SIGNATURES",j,n,Cputime(tt);if n eq 23 then break;end if;
end for;
assert #Seqset(sigs) eq 23;groups:=AssociativeArray();for i in [1..#ps] do if IsDefined(groups,sigs[i]) then Append(~groups[sigs[i]],i);else groups[sigs[i]]:=[i];end if;end for;
blocks:=[groups[s]:s in Sort(Setseq(Keys(groups)))];assert Sort([#v:v in blocks]) eq [1:i in [1..20]] cat [2,2,2];delta_blocks:=[];
for p in ps do sig:=[Valuation(u,p):u in ugs];assert IsDefined(groups,sig);Append(~delta_blocks,groups[sig]);end for;
Write("work/deg10/C232_resolution_20261030/hecke5_cusp_geometry.m",Sprint([*geom[1],[EP(p):p in ps],sigs,blocks,delta_blocks,[EF(u):u in us],[EF(u):u in ugs]*]):Overwrite:=true);print "C232_CUSP_SIGNATURE_GEOMETRY_COMPLETE",blocks,delta_blocks,Cputime(tt);quit;
