SetColumns(0);SetNthreads(1);SetQuitOnError(true);SetMemoryLimit(9*10^9);Alarm(950);tt:=Cputime();root:="work/deg10/continuation_20261031a/C70/";
packets:=[eval Read(root cat name cat ".m"):name in INPUTNAMES];p:=packets[1];assert &and[q[1] eq p[1] and q[2] eq p[2]:q in packets];
k:=GF(3);PP:=PolynomialRing(k);K<t>:=FunctionField(k);R<X>:=PolynomialRing(K);
function DR(v)return K!(PP!v[1])/K!(PP!v[2]);end function;
F<y>:=FunctionField(R![DR(v):v in p[1]]);assert Degree(F) eq 15 and Genus(F) eq 25;
function DF(v)return F![DR(c):c in v];end function;
xy:=[DF(v):v in p[2]];OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);cache:=AssociativeArray();
function ParseP(row)
 O:=row[2] eq 1 select OF else OI;I:=ideal<O|[DF(v):v in row[3]]>;assert IsPrime(I);P:=Place(I);assert Degree(P) eq row[1];return P;
end function;
function Key2(v)
 L:=Parent(v[1]);d:=LCM([Degree(MinimalPolynomial(s,k)):s in v]);found:=false;
 for bb in[[1,0],[0,1],[1,1],[1,2]] do
  a:=&+[L!bb[j]*v[j]:j in[1..2]];mp:=MinimalPolynomial(a,k);
  if Degree(mp) eq d then bits:=bb;found:=true;break;end if;
 end for;assert found;B:=Matrix(k,[Eltseq(a^j,k):j in[0..d-1]]);assert Rank(B) eq d;co:=[Eltseq(Solution(B,Vector(k,Eltseq(s,k)))):s in v];
 assert &and[&+[L!co[j][i+1]*a^i:i in[0..d-1]] eq v[j]:j in[1..2]];
 return [*d,bits,[Integers()!v:v in Coefficients(mp)],[[Integers()!v:v in r]:r in co]*];
end function;
function KP(P)return Key2([Evaluate(v,P):v in xy]);end function;
out:=[* *];nm:=0;np:=0;nc:=0;
for packet in packets do for row in packet[3] do
 ci,d,maps:=Explode(row);rs:=[* *];good:=false;
 for mi in[1..#maps] do
  r:=maps[mi];ps:=ParseP(r[1]);pd:=ParseP(r[3]);np+:=2;assert d eq r[2]*Degree(ps) and 2*d eq r[4]*Degree(pd);T:=DivisorGroup(F)!0;obs:=[* *];assert &+[q[1]:q in r[5]] eq 3;
  for q in r[5] do pt:=ParseP(q[3]);np+:=1;assert d*q[1] eq q[2]*Degree(pt);T+:=q[2]*Divisor(pt);Append(~obs,[*q[1],q[2],KP(pt)*]);end for;
  assert Degree(T) eq 3*d;B:=T-r[2]*Divisor(ps)-r[4]*Divisor(pd);assert Degree(B) eq 0;hh:=#Basis(B:Reduction:=false,Simplification:="None");assert hh eq r[6];if hh eq 0 then good:=true;end if;
  Append(~rs,[*mi,KP(ps),r[2],KP(pd),r[4],obs,hh*]);nm+:=1;
 end for;assert good;Append(~out,[*ci,d,rs*]);nc+:=1;
 if nc mod 10 eq 0 then print "C70_PARENT_NATIVE_PROGRESS",TAG,nc,nm,Cputime(tt);end if;
end for;end for;
Write(root cat TAG cat ".m",Sprint(out):Overwrite:=true);print "C70_PARENT_NATIVE_COMPLETE",TAG,nc,nm,np,Cputime(tt);quit;
