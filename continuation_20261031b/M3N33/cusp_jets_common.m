SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(8*10^9);Alarm(1100);tt:=Cputime();root:="work/deg10/continuation_20261031b/M3N33/";dat:=eval Read(root cat "relative_kummer_frame.m");bps:=eval Read(root cat "relative_base_cusps.m");F2:=GF(2);R2:=PolynomialRing(F2);k<w>:=ext<F2|R2!dat[1]>;
function DC(a)return k!(a mod 2)+k!(a div 2)*w;end function;
function EC(a)v:=Eltseq(k!a,F2);return Integers()!v[1]+2*Integers()!v[2];end function;
B<u>:=FunctionField(k);PB:=PolynomialRing(k);RY<Y>:=PolynomialRing(B);
function DR(v)return B!(PB![DC(c):c in v[1]])/(B!(PB![DC(c):c in v[2]]));end function;
F<v>:=FunctionField(RY![DR(c):c in dat[2]]);RZ<Z>:=PolynomialRing(F);
function DF(a)return F![DR(c):c in a];end function;
L<r>:=FunctionField(RZ![DF(c):c in dat[3]]:Check:=false);assert Genus(L) eq 81;
function DL(a)return L![DF(c):c in a];end function;
OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);OL:=MaximalOrderFinite(L);OIL:=MaximalOrderInfinite(L);
function BasePrime(a)O:=a[2] eq 1 select OF else OI;I:=ideal<O|[DF(c):c in a[3]]>;assert IsPrime(I);P:=Place(I);assert Degree(P) eq a[1];return P;end function;
function ER(a)return [[EC(c):c in Coefficients(Numerator(a))],[EC(c):c in Coefficients(Denominator(a))]];end function;
function EF(a)return [ER(c):c in Eltseq(F!a)];end function;
function EL(a)return [EF(c):c in Eltseq(L!a)];end function;
function EP(P)a,b:=TwoGenerators(P);assert ideal<Order(Ideal(P))|a,b> eq Ideal(P);return [*Degree(P),IsFinite(P) select 1 else 0,[EL(a),EL(b)]*];end function;
bb:=[DL(c):c in dat[5]];assert #bb eq 81;ur:=#dat[8][3] div 81;U:=Matrix(k,ur,81,[DC(c):c in dat[8][3]]);assert Rank(U) eq ur;du:=Differential(L!u);rows:=[* *];basesdone:=0;
print "M333_CUSP_JET_SHARD_MODEL_READY",SHARD,ur,#bps,Cputime(tt);
for ii in [i:i in[1..#bps]|(i-1) mod 2 eq SHARD] do
 P:=BasePrime(bps[ii]);baseunif:=UniformizingElement(P);a,b:=TwoGenerators(P);raw:=Decomposition(L,P);qs:=[Q:Q in raw|Valuation(L!a,Q) gt 0 and Valuation(L!b,Q) gt 0];assert &+[Valuation(L!baseunif,Q)*Degree(Q):Q in qs] eq 3*Degree(P);
 for jj in[1..#qs] do
  Q:=qs[jj];d:=Degree(Q);assert d in[1,5];m:=5 div d;z:=UniformizingElement(Q);assert Valuation(z,Q) eq 1;dz:=Differential(z)/du;functions:=[b/dz:b in bb];assert &and[b eq 0 or Valuation(b,Q) ge 0:b in functions];
  if d eq 5 then K:=ResidueClassField(Q);vals:=[[K!Evaluate(b,Q)]:b in functions];else CP,cp:=Completion(L,Q:Precision:=30);K:=BaseRing(CP);series:=[cp(b):b in functions];vals:=[[Coefficient(s,j):j in[0..m-1]]:s in series];end if;
  Embed(k,K);J:=Matrix(k,81,m*d,&cat[&cat[Eltseq(K!a,k):a in row]:row in vals]);assert Rank(J) eq m*d;UJ:=U*J;Append(~rows,[*ii,jj,EP(Q),EL(z),EL(dz),d,m,[[EC(a):a in row]:row in Rows(J)],[[EC(a):a in row]:row in Rows(UJ)], [EC(a):a in Coefficients(DefiningPolynomial(K,k))]*]);
  Write(root cat Sprintf("cusp_jets_s%o_partial.m",SHARD),Sprint(rows):Overwrite:=true);print "M333_CUSP_JET_READY",SHARD,ii,jj,d,m,Rank(UJ),Cputime(tt);
 end for;basesdone+:=1;
end for;
Write(root cat Sprintf("cusp_jets_s%o.m",SHARD),Sprint(rows):Overwrite:=true);print "D31B_M333_CUSP_JET_SHARD_COMPLETE",SHARD,basesdone,#rows,Cputime(tt);quit;
