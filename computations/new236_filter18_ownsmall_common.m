SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(5*10^9);tt:=Cputime();k:=GF(5);PP:=PolynomialRing(k);
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
cache:=AssociativeArray();images:=[* *];
procedure Classes(~cl,~cache,~images,ky,d,K)
 z:=K.1;b:=Evaluate(PP!ky[1..d],z);c:=Evaluate(PP!ky[d+1..2*d],z);E:=EllipticCurve([K|0,c,0,b*(1-b-c),0]);P:=E![0,0];Q:=E![b,b];assert Order(Q) eq 36 and P ne 18*Q;cl:=[];
 for mapno in [1,2] do
  if mapno eq 1 then EG:=E;U:=P;V:=2*Q;
  else T:=18*Q;EG,iso:=IsogenyFromKernel(E,PolynomialRing(K).1-T[1]);U:=iso(P);V:=iso(Q);end if;assert Order(V) eq 18 and U ne 9*V;bc:=BC(EG,U,V);kk:=[d] cat Key(bc,d);
  if IsDefined(cache,kk) then j:=cache[kk];answer:=q(Z!images[j][4]);
  else
   fast,uv:=Pair(bc);pl:=C48PairAFPlaceLinear(fast select coordsUV else coordsBC,uv);e:=Degree(pl);assert d mod e eq 0;mult:=d div e;v:=phi(Divisor(pl)-e*A);answer:=mult*v;I:=Ideal(pl);flag:=Order(I) eq OF select 1 else 0;
   Append(~images,[*kk,e,mult,Eltseq(answer),flag,[EF(v):v in Basis(I)],Eltseq(v)*]);j:=#images;cache[kk]:=j;
  end if;Append(~cl,Eltseq(answer));if d ge 8 and not answer in H then break;end if;
 end for;
end procedure;
for d in DS do
 src:=eval Read(Sprintf("work/deg10/C236_resolution_20261030/enumeration/source%o.m",d));assert src[1] eq d and src[3];K<zz>:=ext<k|PP!src[2]>;own:=[* *];small:=[* *];allkeys:={};
 oi:=0;for rec in src[9] do
  if rec[1] ne d then continue;end if;oi+:=1;if (oi-1) mod SHARDS ne SHARD then continue;end if;ky:=rec[2];cl:=[];Classes(~cl,~cache,~images,ky,d,K);Append(~own,[*d,oi,ky,cl*]);
  if false then
   b:=Evaluate(PP!ky[1..d],zz);c:=Evaluate(PP!ky[d+1..2*d],zz);E:=EllipticCurve([K|0,c,0,b*(1-b-c),0]);P:=E![0,0];Q:=E![b,b];ims:={Key(BC(E,P+18*bb*Q,cc*P+u*Q),d):bb,cc in [0,1],u in [1,5,7,11,13,17]};assert IsDisjoint(allkeys,ims);allkeys join:=ims;
  end if;
  if #own mod 1000 eq 0 then print "C236_FILTER18_OWN_PROGRESS",TAG,d,#own,#images,Cputime(tt);end if;
 end for;
 if false then
  for ky in Sort(Setseq(allkeys)) do cl:=[];Classes(~cl,~cache,~images,ky,d,K);Append(~small,[*d,#small+1,ky,cl*]);end for;
 end if;
 Write(Sprintf("work/deg10/C236_resolution_20261030/enumeration/filter18_D%o_%o.m",d,TAG),Sprint([*src[2],own,small*]):Overwrite:=true);print "C236_FILTER18_DEGREE_COMPLETE",d,#own,#small,#images,Cputime(tt);
end for;
Write(Sprintf("work/deg10/C236_resolution_20261030/enumeration/filter18_images_%o.m",TAG),Sprint(images):Overwrite:=true);fprintf witness,"]\n";Flush(witness);print "C236_FILTER18_BATCH_COMPLETE",TAG,#images,Cputime(tt);quit;
