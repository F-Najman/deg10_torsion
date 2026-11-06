SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(5*10^9);Alarm(1100);tt:=Cputime();AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");load "work/deg10/continuation_20261031b/C60/inputs/portable_forward.m";load "work/deg10/continuation_20261031b/C60/inputs/pair_af_place_linear.m";load "work/deg10/continuation_20261031b/C60/exceptional_fibres_data.m";
root:="work/deg10/continuation_20261031b/C60/";a:=eval Read(root cat "target30_proven_picard.m");k:=GF(7);PK:=PolynomialRing(k);B<t>:=FunctionField(k);R<Y>:=PolynomialRing(B);
function Rat(v)return B!(PK!v[1])/B!(PK!v[2]);end function;
F<y>:=FunctionField(R![Rat(v):v in a[1]]);assert Genus(F)eq 9;
function FF(v)return F![Rat(z):z in v];end function;
OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
function Prime(v)O:=v[2]eq 1 select OF else OI;I:=ideal<O|[FF(z):z in v[3]]>;p:=Place(I);assert Degree(p)eq v[1];return p;end function;
fb:=[Prime([*a[3][i],a[2][i],a[4][i]*]):i in[1..#a[3]]];A:=Divisor(Prime(a[12]));Z:=FreeAbelianGroup(#a[6][1]);G,q:=quo<Z|[Z!v:v in a[6]]>;coords:=[FF(v):v in a[7]];ph,witness:=C240PortableDirectRR(F,fb,[q(Z!v):v in a[5]],G,A,root cat TAG cat "_class.log");
function Key(values)
 K:=Parent(values[1]);e:=LCM([Degree(MinimalPolynomial(v,k)):v in values]);found:=false;for bits in[[1,0],[0,1]]cat[[1,z]:z in[1..6]]do gen:=&+[K!bits[i]*values[i]:i in[1,2]];pol:=MinimalPolynomial(gen,k);if Degree(pol)eq e then found:=true;chosen:=bits;break;end if;end for;assert found;M:=Matrix(k,[Eltseq(gen^i,k):i in[0..e-1]]);cs:=[Eltseq(Solution(M,Vector(k,Eltseq(v,k)))):v in values];return [e]cat chosen cat[Integers()!z:z in Coefficients(pol)]cat&cat[[Integers()!z:z in row]:row in cs];
end function;
images:=[* *];cache:=AssociativeArray();records:=[* *];for ci in[SHARD+1..#C230Points by NSHARDS]do d,bits,pol,co:=Explode(C230Points[ci]);K<z>:=ext<k|PK!pol>;bc:=[Evaluate(PK!c,z):c in co];b,c:=Explode(bc);E:=EllipticCurve([K|0,c,0,b*(1-b-c),0]);P2:=E![0,0];Q:=E![b,b];assert Order(P2)eq 2 and Order(Q)eq 30 and P2 ne 15*Q;qs:=[Q,P2+16*Q,P2+Q];rows:=[* *];
 for v in[1..3]do assert Order(qs[v])eq 30;xy:=TateNormalForm_xy(E,qs[v]);key:=Key(xy);e:=key[1];assert d mod e eq 0;if IsDefined(cache,key)then ix:=cache[key];cl0:=q(Z!images[ix][2]);else pt:=C48PairAFPlaceLinear(coords,xy);assert Degree(pt)eq e;cl0:=ph(Divisor(pt)-e*A);Append(~images,[*key,Eltseq(cl0),e*]);ix:=#images;cache[key]:=ix;end if;cl:=(d div e)*cl0;Append(~rows,[*v,ix,Eltseq(cl),e,d div e*]);end for;Append(~records,[*ci,d,rows*]);if #records mod 50 eq 0 then Write(root cat TAG cat "_progress.m",Sprint([*records,images*]):Overwrite:=true);print "D31B_C60_THREE_MAPS_PROGRESS",SHARD,#records,#images,Cputime(tt);end if;
end for;
Write(root cat TAG cat ".m",Sprint([*records,images*]):Overwrite:=true);fprintf witness,"]\n";Flush(witness);print "D31B_C60_THREE_MAPS_COMPLETE",SHARD,#records,#images,Cputime(tt);quit;
