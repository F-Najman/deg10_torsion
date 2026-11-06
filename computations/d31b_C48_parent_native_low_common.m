SetColumns(0);SetNthreads(1);SetQuitOnError(true);SetMemoryLimit(8*10^9);Alarm(1100);tt:=Cputime();root:="work/deg10/continuation_20261031b/C48/";k:=GF(5);PK:=PolynomialRing(k);B<t>:=FunctionField(k);R<Y>:=PolynomialRing(B);md:=eval Read(root cat Sprintf("hecke%o_s0_model.m",Q));
function DR(a)return B!(PK!a[1])/(B!(PK!a[2]));end function;
f:=R![DR(a):a in md[1]];assert IsIrreducible(f);F<u>:=FunctionField(f);assert Genus(F)eq 37;OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
function DF(v)return F![DR(a):a in v];end function;xy:=[DF(v):v in md[2]];
function Prime(pd)O:=pd[2]eq 1 select OF else OI;gs:=[DF(v):v in pd[3]];I:=ideal<O|gs>;assert IsIntegral(I)and IsPrime(I);P:=Place(I);assert Degree(P)eq pd[1]and Ideal(P)eq I;return P;end function;
function Key(P)
 v:=[Evaluate(z,P):z in xy];L:=Parent(v[1]);dd:=Degree(P);assert LCM([Degree(MinimalPolynomial(z,k)):z in v])eq dd;found:=false;
 for bits in[[1,0],[0,1],[1,1],[1,2],[1,3],[1,4]]do g:=&+[L!bits[j]*v[j]:j in[1,2]];pol:=MinimalPolynomial(g,k);if Degree(pol)eq dd then found:=true;chosen:=bits;break;end if;end for;assert found;M:=Matrix(k,[Eltseq(g^i,k):i in[0..dd-1]]);assert Rank(M)eq dd;co:=[Eltseq(Solution(M,Vector(k,Eltseq(z,k)))):z in v];return[*dd,chosen,[Integers()!a:a in Coefficients(pol)],[[Integers()!a:a in r]:r in co]*];
end function;
packet:=eval Read(root cat Sprintf("low_hecke%o_merged.m",Q));points:=packet[3];cases:=packet[4];cache:=AssociativeArray();out:=[* *];np:=0;
for row in points do id,d,q,ps,ns,pd,nd,hec:=Explode(row);assert q eq Q;P:=Prime(ps);V:=Prime(pd);np+:=2;assert ns*Degree(P)eq d and nd*Degree(V)eq Q*d;A:=-ns*Divisor(P)-nd*Divisor(V);hh:=[* *];assert &+[h[1]:h in hec]eq Q+1;
 for h in hec do T:=Prime(h[3]);np+:=1;assert h[2]*Degree(T)eq d*h[1];A+:=h[2]*Divisor(T);Append(~hh,[*h[1],h[2],Key(T)*]);end for;assert Degree(A)eq 0;cache[id]:=A;Append(~out,[*id,d,Q,Key(P),ns,Key(V),nd,hh*]);if #out mod 10 eq 0 then print "C48_PARENT_NATIVE_LOW_POINT",Q,#out,np,Cputime(tt);end if;
end for;
rr:=[* *];for cas in cases do id,parts,ids,cusp,oldh:=Explode(cas);assert &+parts+cusp eq 10 and #parts eq #ids;D:=&+[cache[id]:id in ids];assert Degree(D)eq 0;h0:=#Basis(D:Reduction:=false,Simplification:="None");assert h0 eq 0 and oldh eq 0;Append(~rr,[*id,parts,ids,cusp,h0*]);if #rr mod 10 eq 0 then print "C48_PARENT_NATIVE_LOW_CASE",Q,#rr,Cputime(tt);end if;end for;
Write(root cat Sprintf("parent_native_low_q%o.m",Q),Sprint([*out,rr*]):Overwrite:=true);print "D31B_C48_PARENT_NATIVE_LOW_COMPLETE",Q,#out,#rr,np,Cputime(tt);quit;
