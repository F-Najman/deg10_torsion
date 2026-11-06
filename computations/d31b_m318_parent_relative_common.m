SetColumns(0);SetNthreads(1);SetQuitOnError(true);SetMemoryLimit(8*10^9);Alarm(1000);tt:=Cputime();root:="work/deg10/continuation_20261031b/M3N18/";name:=PRIME eq 7 select "relative_cubic_rank.m"else"relative_cubic_rank_p13.m";a:=eval Read(root cat name);assert #a[6]eq 0;k:=GF(PRIME);PK:=PolynomialRing(k);B<u>:=FunctionField(k);R<Y>:=PolynomialRing(B);
function DR(v)return B!(PK!v[1])/(B!(PK!v[2]));end function;
f:=R![DR(v):v in a[1]];assert IsIrreducible(f);F<v>:=FunctionField(f);assert Genus(F)eq 10;OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
function DF(vv)return F![DR(z):z in vv];end function;
function Prime(pd)O:=pd[2]eq 1 select OF else OI;I:=ideal<O|[DF(z):z in pd[3]]>;assert IsPrime(I);P:=Place(I);assert Degree(P)eq pd[1];return P;end function;
cs:=[Prime(p):p in a[3]];du:=Differential(F!u);bf:=[DF(z):z in a[2]];HH,hm:=SpaceOfHolomorphicDifferentials(F);assert Dimension(HH)eq 10;M:=Matrix(k,[Eltseq((z*du)@@hm):z in bf]);assert Rank(M)eq 10;
function Coords(z)return Eltseq(Solution(M,Vector(Eltseq(z@@hm))));end function;
out:=[* *];nc:=0;
for block in a[5]do
 ii,par,pn,imgs,oldTM,rows:=Explode(block);N:=Prime(pn);assert Degree(N)eq 3;up:=DF(imgs[2]);vp:=DF(imgs[1]);bm:=hom<B->F|up>;tau:=hom<F->F|bm,vp>;assert tau(tau(tau(F.1)))eq F.1 and tau(tau(tau(F!u)))eq F!u and(tau(F.1)ne F.1 or tau(F!u)ne F!u);gg:=pn[3];assert &and[Valuation(tau(DF(z)),N)gt 0:z in gg];dt:=Differential(tau(F!u))/du;TM:=Matrix(k,[Coords(tau(z)*dt*du):z in bf]);assert TM eq Matrix(k,oldTM);U:=RowSpace(TM-1);assert Dimension(U)eq 6;
 #rows eq 171;rr:=[* *];
 for r in rows do i,j,rk,oldW:=Explode(r);D:=Divisor(N)+Divisor(cs[i])+Divisor(cs[j]);assert Degree(D)eq 5;V,vm:=DifferentialSpace(D);W:=sub<VectorSpace(k,10)|[Coords(vm(z)):z in Basis(V)]>;assert W eq sub<VectorSpace(k,10)|[Vector(k,z):z in oldW]>;rank:=Dimension(U)-Dimension(U meet W);assert rank eq rk and rank eq 4;Append(~rr,[i,j,rank,Dimension(W)]);nc+:=1;end for;
 Append(~out,[*ii,par,rr*]);print "M318_PARENT_RELATIVE_BLOCK",PRIME,ii,nc,Cputime(tt);
end for;
assert nc eq(PRIME eq 7 select 684 else 1710);Write(root cat Sprintf("parent_relative_cubic_p%o.m",PRIME),Sprint(out):Overwrite:=true);print "D31B_M318_PARENT_RELATIVE_COMPLETE",PRIME,#out,nc,Cputime(tt);quit;
