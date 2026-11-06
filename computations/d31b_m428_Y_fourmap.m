SetColumns(0);Alarm(1100);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(6*10^9);tt:=Cputime();AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");load "work/deg10/continuation_20261031b/M4N28/inputs/portable_forward.m";load "work/deg10/continuation_20261031b/M4N28/inputs/pair_af_place_linear.m";
root:="work/deg10/continuation_20261031b/M4N28/";a:=eval Read(root cat "inputs/target28_proven_picard.m");k:=GF(3);PK:=PolynomialRing(k);B<t>:=FunctionField(k);R<Y>:=PolynomialRing(B);
function Rat(v)return B!(PK!v[1])/B!(PK!v[2]);end function;
F<y>:=FunctionField(R![Rat(v):v in a[1]]);assert Genus(F)eq 10;
function FF(v)return F![Rat(z):z in v];end function;
OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
function Prime(v)O:=v[2]eq 1 select OF else OI;I:=ideal<O|[FF(z):z in v[3]]>;p:=Place(I);assert Degree(p)eq v[1];return p;end function;
fb:=[Prime([*a[3][i],a[2][i],a[4][i]*]):i in[1..#a[3]]];A:=Divisor(Prime(a[12]));Z:=FreeAbelianGroup(#a[6][1]);G,q:=quo<Z|[Z!v:v in a[6]]>;fc:=[q(Z!v):v in a[5]];H:=sub<G|[q(Z!v):v in a[13]]>;GG,pi:=quo<G|H>;coords:=[FF(v):v in a[7]];ph,witness:=C240PortableDirectRR(F,fb,fc,G,A,root cat "fourmap_native_class.log");
allowed:=AssociativeArray();for d in[0..5]do allowed[d]:={GG|};end for;allowed[0]:={GG!0};for i in[1..#a[9]]do w:=a[9][i];for d in[w..5]do allowed[d]join:={z+pi(q(Z!a[10][i])):z in allowed[d-w]};end for;end for;
print "D31B_M428_FOURMAP_GROUP_READY",Invariants(G),Invariants(H),[#allowed[d]:d in[0..5]],Cputime(tt);

Attach("work/deg10/continuation_20261031b/M4N28/inputs/fast_diamond_hash_m2.m");Attach("work/deg10/continuation_20261031b/M4N28/inputs/fast_full_aut_hash_m2_28.m");
X:=MDX11(2,28,k);assert Genus(X)eq 25;tmp:=eval Read("work/deg10/continuation_20261031b/C56/projected_Y28.m");fd:=tmp[1];fields:=AssociativeArray();for z in fd do fields[z[1]]:=ext<k|PK!z[2]>;end for;
cache:=AssociativeArray();images:=[* *];records:=[* *];left:=[* *];coverage:=[* *];
for d in[5,6,8,10]do
 xs:=NoncuspidalPlaces(X,d);rs,labels,allkeys:=PlacesUpToFullAutM2TwentyEight(X,xs);K:=fields[d];ct:=0;
 for ri in[1..#rs]do
  bc:=MDM2ToCommonField([Evaluate(v,rs[ri]):v in X`_coordinates],K);b,c:=Explode(bc);E:=EllipticCurve([K|0,c,0,b*(1-b-c),0]);P2:=E![0,0];Q:=E![b,b];assert Order(Q)eq 28 and 2*P2 eq E!0 and P2 ne E!0 and P2 ne 14*Q;rr:=[v[1]:v in Roots(DivisionPolynomial(E,2))];assert #rr eq 3;
  full4:=IsEven(d)and &and[IsSquare(rr[i]-rr[j]):i,j in[1..3]|i lt j];if d ne 5 and not full4 then continue;end if;ct+:=1;mul:=d eq 5 select 2 else 1;degpart:=d*mul;ky:=MDM2FrobeniusKey(bc,K,3,d);cls:=[* *];obs:=[* *];
  for v in[1..4]do
   if v le 2 then ET:=E;QT:=v eq 1 select Q else P2+Q;else S:=v eq 3 select P2 else P2+14*Q;ET,ff:=IsogenyFromKernel(E,PolynomialRing(K).1-S[1]);QT:=ff(Q);end if;assert Order(QT)eq 28;xy:=TateNormalForm_xy(ET,QT);ck:=[d]cat MDM2FrobeniusKey(xy,K,3,d);
   if IsDefined(cache,ck)then ix:=cache[ck];cl:=q(Z!images[ix][2]);e:=images[ix][3];else pt:=C48PairAFPlaceLinear(coords,xy);e:=Degree(pt);assert d mod e eq 0;cl:=(d div e)*ph(Divisor(pt)-e*A);Append(~images,[*ck,Eltseq(cl),e*]);ix:=#images;cache[ck]:=ix;end if;Append(~cls,Eltseq(mul*cl));Append(~obs,[v,ix,e,mul*d div e]);
  end for;flags:=[-pi(q(Z!cl))in allowed[10-degpart]:cl in cls];si:=#records+1;Append(~records,[*si,d,mul,ky,cls,obs,flags*]);if &and flags then Append(~left,[*si,[si],10-degpart*]);end if;print "D31B_M428_Y_FOURMAP_CASE",d,ri,si,flags,Cputime(tt);
 end for;Append(~coverage,[*d,#xs,#rs,ct*]);print "D31B_M428_Y_FOURMAP_DEGREE",d,#xs,#rs,ct,Cputime(tt);
end for;
fprintf witness,"]\n";Flush(witness);Write(root cat "fourmap_Y28.m",Sprint([*fd,records,left,coverage,images,[Eltseq(v):v in Rows(RelationMatrix(G))],a[13]*]):Overwrite:=true);print "D31B_M428_Y_FOURMAP_COMPLETE",#records,#left,coverage,Cputime(tt);quit;
