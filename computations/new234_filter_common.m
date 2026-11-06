SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(4*10^9);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
load "work/deg10/C234_resolution_20261030/deps/pair_af_place_linear.m";
load "work/deg10/C234_resolution_20261030/deps/portable_forward.m";
a:=eval Read("work/deg10/C234_resolution_20261030/target_proven_picard.m");k:=GF(3);PP:=PolynomialRing(k);B<t>:=FunctionField(k);RF<Y>:=PolynomialRing(B);
function Rat(v) return B!(PP!v[1])/B!(PP!v[2]);end function;
f:=RF![Rat(v):v in a[1]];f/:=LeadingCoefficient(f);F<y>:=FunctionField(f);assert Genus(F) eq 21;
function GetFunction(data) return F![Rat(v):v in data];end function;
OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);fb:=[];
for j in [1..#a[4]] do OO:=a[2][j] eq 1 select OF else OI;I:=ideal<OO|[GetFunction(v):v in a[4][j]]>;pl:=Place(I);assert Degree(pl) eq a[3][j];Append(~fb,pl);end for;
Z:=FreeAbelianGroup(#a[6][1]);G,q:=quo<Z|[Z!v:v in a[6]]>;fc:=[q(Z!v):v in a[5]];A:=Divisor(fb[a[8]]);coords:=[GetFunction(v):v in a[7]];
phi,witness:=C240PortableDirectRR(F,fb,fc,G,A,Sprintf("work/deg10/C234_resolution_20261030/filter%o_RR.log",TAG));
hg:=[q(Z!v):v in a[9]];H:=sub<G|hg>;assert Invariants(H) eq [8760,595680];QH,qh:=quo<G|H>;assert #TorsionSubgroup(QH) eq 128;
assert &and[q(Z!v) in H:v in a[11]];print "C234_FILTER_READY",TAG,Invariants(G),Cputime(tt);
function Key(v,D)
 keys:=[];for i in [1..D] do Append(~keys,&cat[[Integers()!a:a in Eltseq(z,k)]:z in v]);v:=[z^3:z in v];end for;return Sort(keys)[1];end function;
function BC(E,P,Q) ai:=aInvariants(E);xx:=Q[1]-P[1];yy:=Q[2];return [xx^3/yy^2,(ai[2]+3*P[1])*xx^2/yy^2];end function;
function ER(v) return [[Integers()!c:c in Coefficients(Numerator(v))],[Integers()!c:c in Coefficients(Denominator(v))]];end function;
function EF(v) return [ER(c):c in Eltseq(F!v)];end function;
fields:=AssociativeArray();own:=AssociativeArray();small:=AssociativeArray();source:=[* *];coverage:=[* *];
for d in [4..10] do
 src:=eval Read(Sprintf("work/deg10/C234_resolution_20261030/source%o.m",d));assert src[1] eq d and src[3];K<z>:=ext<k|PP!src[2]>;fields[d]:=K;own[d]:=[* *];keys:={};
 for r in src[7] do
  if r[1] ne d then continue;end if;row:=[*d,#own[d]+1,r[2]*];Append(~own[d],row);Append(~source,row);
  if d le 6 then
   ky:=r[2];b:=Evaluate(PP!ky[1..d],z);c:=Evaluate(PP!ky[d+1..2*d],z);E:=EllipticCurve([K|0,c,0,b*(1-b-c),0]);P:=E![0,0];Q:=E![b,b];ts:=[P,17*Q,P+17*Q];V:=18*Q;ims:={Key(BC(E,ts[i],ts[j]+a*V),d):i,j in [1..3],a in [1..8]|i ne j};assert IsDisjoint(keys,ims);keys join:=ims;
  end if;
 end for;
 if d le 6 then small[d]:=[*[*d,i,ky*]:i->ky in Sort(Setseq(keys))*];assert #small[d] eq src[6][#Divisors(d)];end if;
 Append(~coverage,[*d,src[2],#own[d]*]);
end for;
supports:=[*[1,r[1],r[2]]:r in source*];for ty in [[4,4],[4,5],[4,6],[5,5]] do for i in [1..#own[ty[1]]],j in [1..#small[ty[2]]] do Append(~supports,[2,ty[1],i,ty[2],j]);end for;end for;
cache:=AssociativeArray();images:=[* *];tests:=[* *];remaining:=[* *];
procedure Push(~answer,~cache,~images,row,variant)
 d:=row[1];ky:=row[3];K:=fields[d];z:=K.1;b:=Evaluate(PP!ky[1..d],z);c:=Evaluate(PP!ky[d+1..2*d],z);E:=EllipticCurve([K|0,c,0,b*(1-b-c),0]);P:=E![0,0];Q:=E![b,b];assert Order(Q) eq 34 and P ne 17*Q;
 if variant le 2 then EG:=E;U:=variant eq 1 select Q else P+Q;else T:=variant eq 3 select P else P+17*Q;EG,iso:=IsogenyFromKernel(E,PolynomialRing(K).1-T[1]);U:=iso(Q);end if;assert Order(U) eq 34;xy:=TateNormalForm_xy(EG,U);key:=[d] cat Key(xy,d);
 if IsDefined(cache,key) then ti:=cache[key];answer:=q(Z!images[ti][4]);return;end if;
 pl:=C48PairAFPlaceLinear(coords,xy);e:=Degree(pl);assert d mod e eq 0;v:=phi(Divisor(pl)-e*A);answer:=(d div e)*v;I:=Ideal(pl);flag:=Order(I) eq OF select 1 else 0;Append(~images,[*key,e,d div e,Eltseq(answer),flag,[EF(x):x in Basis(I)],Eltseq(v)*]);cache[key]:=#images;
end procedure;
for j in [1..#supports] do
 if (j-1) mod SHARDS ne SHARD then continue;end if;desc:=supports[j];sp:=[*own[desc[2]][desc[3]]*];if desc[1] eq 2 then Append(~sp,small[desc[4]][desc[5]]);end if;cd:=10-&+[row[1]:row in sp];assert cd in [0..6];checks:=[* *];passed:=true;
 for variant in [1..4] do cl:=G!0;for row in sp do v:=G!0;Push(~v,~cache,~images,row,variant);cl+:=v;end for;yes:=cl in H;Append(~checks,[*variant,Eltseq(cl),yes*]);if not yes then passed:=false;break;end if;end for;
 Append(~tests,[*j,desc,cd,checks,passed*]);if passed then Append(~remaining,desc);end if;
 if #tests mod 20 eq 0 then Write(Sprintf("work/deg10/C234_resolution_20261030/filter%o.m",TAG),Sprint([*source,[*[*d,small[d]*]:d in [4..6]*],coverage,images,tests,remaining*]):Overwrite:=true);print "C234_FILTER_PROGRESS",TAG,j,#remaining,Cputime(tt);end if;
end for;
Write(Sprintf("work/deg10/C234_resolution_20261030/filter%o.m",TAG),Sprint([*source,[*[*d,small[d]*]:d in [4..6]*],coverage,images,tests,remaining*]):Overwrite:=true);fprintf witness,"]\n";Flush(witness);print "C234_FILTER_COMPLETE",TAG,#supports,#tests,#images,#remaining,Cputime(tt);quit;
