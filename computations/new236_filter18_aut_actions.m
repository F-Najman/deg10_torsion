TAG:="autactions";SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(5*10^9);tt:=Cputime();k:=GF(5);PP:=PolynomialRing(k);
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
Z0:=FreeAbelianGroup(Ngens(G)-1);G0,p0:=quo<Z0|[Z0!(r[1..Ngens(G)-1]):r in a[6]]>;H0:=sub<G0|[p0(Z0!v[1..Ngens(G)-1]):v in a[9]]>;Q0,q0:=quo<G0|H0>;assert #Q0 eq 12;
function QClass(v) e:=Eltseq(v);assert e[#e] eq 0;return q0(p0(Z0!e[1..#e-1]));end function;
b,c:=Explode(coordsBC);EE:=EllipticCurve([F|0,c,0,b*(1-b-c),0]);jj:=jInvariant(EE);baseclass:=phi(A);places:=[p:p in fb|Valuation(jj,p) ge 0];old:=[];new:=[*[],[],[]*];tests:=[* *];
for pl in places do
 vals:=[Evaluate(v,pl):v in coordsBC];K:=Parent(vals[1]);b,c:=Explode(vals);E:=EllipticCurve([K|0,c,0,b*(1-b-c),0]);P:=E![0,0];Q:=E![b,b];assert Order(Q) eq 18 and P ne 9*Q;d:=Degree(pl);cl:=QClass(phi(Divisor(pl)-d*A));Append(~old,cl);row:=[*Index(fb,pl),Eltseq(cl)*];
 for no in [1..3] do
  PP:=no eq 1 select P+9*Q else P;QQ:=no eq 2 select 5*Q else no eq 3 select P+Q else Q;bc:=BC(E,PP,QQ);fast,uv:=Pair(bc);im:=C48PairAFPlaceLinear(fast select coordsUV else coordsBC,uv);assert Degree(im) eq d;v:=QClass(phi(Divisor(im)-d*A));Append(~new[no],v);Append(~row,Eltseq(v));
 end for;Append(~tests,row);
end for;
free:=FreeAbelianGroup(#old);om:=hom<free->Q0|old>;assert Image(om) eq Q0;maps:=[* *];
for no in [1..3] do nm:=hom<free->Q0|new[no]>;assert Kernel(om) eq Kernel(nm);theta:=hom<Q0->Q0|[nm(Q0.i@@om):i in [1..Ngens(Q0)]]>;assert #Kernel(theta) eq 1;Append(~maps,[Eltseq(theta(Q0.i)):i in [1..Ngens(Q0)]]);end for;
out:=[*Invariants(Q0),[Eltseq(r):r in Rows(RelationMatrix(Q0))],[Eltseq(q0(p0(Z0.i))):i in [1..Ngens(Z0)]],maps,tests*];Write("work/deg10/C236_resolution_20261030/enumeration/target18_aut_actions.m",Sprint(out):Overwrite:=true);fprintf witness,"]\n";Flush(witness);print "C236_TARGET18_AUT_ACTIONS_COMPLETE",#places,Invariants(Q0),maps,Cputime(tt);quit;
