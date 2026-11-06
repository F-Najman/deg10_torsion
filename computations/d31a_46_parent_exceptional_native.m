SetColumns(0);SetNthreads(1);SetQuitOnError(true);SetMemoryLimit(8*10^9);Alarm(1600);tt:=Cputime();root:="work/deg10/continuation_20261031a/C46/";
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
md:=eval Read(root cat "parent_model46.m");raw:=eval Read(root cat "parent_exceptional_hecke.m");k:=GF(3);PP:=PolynomialRing(k);K<t>:=FunctionField(k);R<Z>:=PolynomialRing(K);
function DR(v)return K!(PP!v[1])/K!(PP!v[2]);end function;
ff:=R![DR(v):v in md[1]];assert IsIrreducible(ff);F<u>:=FunctionField(ff);assert Genus(F)eq 45;OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
function DF(v)return F![DR(w):w in v];end function;xy:=[DF(v):v in md[2]];
function Prime(pd)O:=pd[2]eq 1 select OF else OI;gs:=[DF(v):v in pd[3]];assert #gs eq 2;I:=ideal<O|gs>;assert IsIntegral(I)and IsPrime(I);P:=Place(I);assert Degree(P)eq pd[1]and Ideal(P)eq I;return P;end function;
function Key2(vals)
 L:=Parent(vals[1]);d:=LCM([Degree(MinimalPolynomial(v,k)):v in vals]);found:=false;
 for bits in[[1,0],[0,1],[1,1],[1,2]]do g:=&+[L!bits[i]*vals[i]:i in[1..2]];mp:=MinimalPolynomial(g,k);if Degree(mp)eq d then chosen:=bits;found:=true;break;end if;end for;assert found;
 M:=Matrix(k,[Eltseq(g^i,k):i in[0..d-1]]);assert Rank(M)eq d;co:=[Eltseq(Solution(M,Vector(k,Eltseq(v,k)))):v in vals];return[*d,chosen,[Integers()!q:q in Coefficients(mp)],[[Integers()!q:q in row]:row in co]*];
end function;
out:=[* *];np:=0;
for caseRow in raw do
 typ,idx,attempts:=Explode(caseRow);ar:=[* *];
 for at in attempts do
  q,comps,oldh0:=Explode(at);D:=DivisorGroup(F)!0;cr:=[* *];
  for row in comps do
   d,mul,pol,aa,qq,pp,dd,hec:=Explode(row);L<z>:=ext<k|PP!pol>;a:=[Evaluate(PP!v,z):v in aa];E:=EllipticCurve(a);Q:=E![Evaluate(PP!v,z):v in qq];assert Order(Q)eq 46;
   P:=Prime(pp);V:=Prime(dd);np+:=2;assert Degree(P)eq d and Degree(V)eq d;kp:=Key2([Evaluate(v,P):v in xy]);kd:=Key2([Evaluate(v,V):v in xy]);assert kp eq Key2(TateNormalForm_xy(E,Q))and kd eq Key2(TateNormalForm_xy(E,q*Q));D-:=mul*(Divisor(P)+q*Divisor(V));tr:=[* *];assert &+[h[1]:h in hec]eq q+1;
   for h in hec do e,n,pd:=Explode(h);T:=Prime(pd);np+:=1;assert d*e eq n*Degree(T);D+:=mul*n*Divisor(T);Append(~tr,[*e,n,Key2([Evaluate(v,T):v in xy])*]);end for;
   Append(~cr,[*d,mul,kp,kd,tr*]);
  end for;
  assert Degree(D)eq 0;h0:=#Basis(D:Reduction:=false,Simplification:="None");assert h0 eq oldh0;Append(~ar,[*q,cr,h0*]);print "C46_PARENT_NATIVE_ATTEMPT",typ,idx,q,h0,Cputime(tt);
 end for;
 assert ar[#ar][3]eq 0;Append(~out,[*typ,idx,ar*]);
end for;
Write(root cat "parent_exceptional_native.m",Sprint(out):Overwrite:=true);print "C46_PARENT_EXCEPTIONAL_NATIVE_COMPLETE",#out,np,Cputime(tt);quit;
