SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(4*10^9);tt:=Cputime();AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
root:="work/deg10/continuation_20261030c/C105/";data:=eval Read(root cat "source105.m");k:=GF(2);PK:=PolynomialRing(k);K<z>:=ext<k|PK!([r[2]:r in data[1]|r[1] eq 10][1])>;rows:=[* *];
function ABS(v,K,M) return [Integers()!a:a in Eltseq(Solution(M,Vector(GF(2),Eltseq(K!v,GF(2)))))];end function;
for i in [1..#data[3]] do
 d,ky:=Explode(data[3][i]);assert d eq 10;b:=Evaluate(PK!ky[1..d],z);c:=Evaluate(PK!ky[d+1..2*d],z);E:=EllipticCurve([K|1-c,-b,-b,0,0]);P:=E![0,0];Q:=3*P;assert Order(P) eq 105 and Order(Q) eq 35;psi:=DivisionPolynomial(E,3);psi/:=LeadingCoefficient(psi);fac:=Factorization(psi);assert &and[v[2] eq 1:v in fac] and &+[Degree(v[1]):v in fac] eq 4;
 for ff in fac do
  pp:=ff[1];e:=Degree(pp);if e eq 1 then L:=K;s:=-Coefficient(pp,0)/Coefficient(pp,1);else L:=ext<K|pp>;s:=L.1;end if;assert Evaluate(pp,s) eq 0;PL<X>:=PolynomialRing(L);EL:=EllipticCurve([L!v:v in aInvariants(E)]);QL:=EL![L!v:v in Eltseq(Q)];PG:=EL![L!v:v in Eltseq(P)];EG,phi:=IsogenyFromKernel(EL,X-s);QT:=phi(QL);assert Order(QT) eq 35;xy:=TateNormalForm_xy(EG,QT);td:=LCM([Degree(MinimalPolynomial(v,k)):v in xy]);assert (10*e) mod td eq 0;n:=10*e div td;
  ag:=PrimitiveElement(L);po:=MinimalPolynomial(ag,k);M:=Matrix(k,[Eltseq(ag^j,k):j in [0..Degree(po)-1]]);assert Rank(M) eq Degree(po);
  Append(~rows,[*i,e,[Integers()!a:a in Coefficients(po)],ABS(L!K.1,L,M),[ABS(v,L,M):v in aInvariants(EL)],[ABS(PG[j],L,M):j in [1,2]],[ABS(QL[j],L,M):j in [1,2]],ABS(s,L,M),[ABS(v,L,M):v in aInvariants(EG)],[ABS(QT[j],L,M):j in [1,2]],[ABS(v,L,M):v in xy],td,n*]);
 end for;
end for;
Write(root cat "kernel3_explicit.m",Sprint(rows):Overwrite:=true);print "C105_EXPLICIT_KERNEL_EXPORT_COMPLETE",#rows,Cputime(tt);quit;
