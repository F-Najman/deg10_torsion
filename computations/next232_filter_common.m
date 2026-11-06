SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(4*10^9);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
load "work/deg10/C232_resolution_20261030/deps/pair_af_place_linear.m";
load "work/deg10/C232_resolution_20261030/deps/portable_forward.m";
a:=eval Read("work/deg10/C232_resolution_20261030/target_proven_picard.m");k:=GF(3);PP:=PolynomialRing(k);B<t>:=FunctionField(k);RF<Y>:=PolynomialRing(B);
function Rat(v) return B!(PP!v[1])/B!(PP!v[2]);end function;
f:=RF![Rat(v):v in a[1]];f/:=LeadingCoefficient(f);F<y>:=FunctionField(f);assert Genus(F) eq 17;
function GetFunction(data) return F![Rat(v):v in data];end function;
OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);fb:=[];
for j in [1..#a[4]] do OO:=a[2][j] eq 1 select OF else OI;I:=ideal<OO|[GetFunction(v):v in a[4][j]]>;pl:=Place(I);assert Degree(pl) eq a[3][j];Append(~fb,pl);end for;
Z:=FreeAbelianGroup(#a[6][1]);G,q:=quo<Z|[Z!v:v in a[6]]>;fc:=[q(Z!v):v in a[5]];A:=Divisor(fb[a[8]]);coords:=[GetFunction(v):v in a[7]];
phi,witness:=C240PortableDirectRR(F,fb,fc,G,A,Sprintf("work/deg10/C232_resolution_20261030/filter%o_RR.log",TAG));
hg:=[q(Z!v):v in a[9]];H:=sub<G|hg>;assert Invariants(H) eq [2,2,2,4,120,11640];QH,qh:=quo<G|H>;assert #TorsionSubgroup(QH) eq 68;
cc:=[q(Z!v):v in a[11]];allowed:=[*{QH!0}*] cat [*{QH|}:i in [1..4]*];for j in [1..#cc] do w:=a[10][j];if w le 4 then for d in [w..4] do allowed[d+1] join:={z+qh(cc[j]):z in allowed[d-w+1]};end for;end if;end for;assert [#co:co in allowed] eq [1,1,3,3,4];
print "C232_FILTER_READY",TAG,Invariants(G),Cputime(tt);
function Key(v,D)
 keys:=[];for i in [1..D] do Append(~keys,&cat[[Integers()!a:a in Eltseq(z,k)]:z in v]);v:=[z^3:z in v];end for;return Sort(keys)[1];
end function;
function ER(v) return [[Integers()!c:c in Coefficients(Numerator(v))],[Integers()!c:c in Coefficients(Denominator(v))]];end function;
function EF(v) return [ER(c):c in Eltseq(F!v)];end function;
cache:=AssociativeArray();images:=[* *];tests:=[* *];remaining:=[* *];
for d in DS do
 s:=eval Read(Sprintf("work/deg10/C232_resolution_20261030/source%o.m",d));assert s[1] eq d and s[3];K<z>:=ext<k|PP!s[2]>;
 for i in [1..#s[7]] do
  if (i-1) mod SHARDS ne SHARD then continue;end if;r:=s[7][i];assert r[1] eq d;ky:=r[2];b:=Evaluate(PP!ky[1..d],z);c:=Evaluate(PP!ky[d+1..2*d],z);E:=EllipticCurve([K|0,c,0,b*(1-b-c),0]);P:=E![0,0];Q:=E![b,b];assert Order(Q) eq 32 and P ne 16*Q;
  checks:=[* *];passed:=true;
  for variant in [1..4] do
   if variant le 2 then EG:=E;U:=variant eq 1 select Q else P+Q;
   else T:=variant eq 3 select P else P+16*Q;EG,iso:=IsogenyFromKernel(E,PolynomialRing(K).1-T[1]);U:=iso(Q);end if;assert Order(U) eq 32;
   xy:=TateNormalForm_xy(EG,U);key:=[d] cat Key(xy,d);
   if IsDefined(cache,key) then ti:=cache[key];e:=images[ti][2];cl:=q(Z!images[ti][4]);
   else
    pl:=C48PairAFPlaceLinear(coords,xy);e:=Degree(pl);assert d mod e eq 0;v:=phi(Divisor(pl)-e*A);cl:=(d div e)*v;
    I:=Ideal(pl);flag:=Order(I) eq OF select 1 else 0;Append(~images,[*key,e,d div e,Eltseq(cl),flag,[EF(x):x in Basis(I)],Eltseq(v)*]);ti:=#images;cache[key]:=ti;
   end if;
   yes:=-qh(cl) in allowed[11-d];Append(~checks,[*variant,ti,e,d div e,Eltseq(cl),yes*]);if not yes then passed:=false;break;end if;
  end for;
  Append(~tests,[*d,i,ky,checks,passed*]);if passed then Append(~remaining,[*d,i,ky*]);end if;
  Write(Sprintf("work/deg10/C232_resolution_20261030/filter%o.m",TAG),Sprint([*images,tests,remaining*]):Overwrite:=true);
  print "C232_FILTER_CASE",TAG,d,i,passed,Cputime(tt);
 end for;
end for;
fprintf witness,"]\n";Flush(witness);print "C232_FILTER_COMPLETE",TAG,#tests,#images,#remaining,Cputime(tt);quit;
