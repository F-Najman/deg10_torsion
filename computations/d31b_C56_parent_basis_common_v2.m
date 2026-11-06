SetColumns(0);SetNthreads(1);SetQuitOnError(true);SetMemoryLimit(8*10^9);Alarm(850);tt:=Cputime();root:="work/deg10/continuation_20261031b/C56/kummer_rr/";cases:=[271,327,491,578,65,818,837,978,1009];k:=GF(3);PK:=PolynomialRing(k);B<t>:=FunctionField(k);R<Y>:=PolynomialRing(B);out:=[* *];
for ci in[1..#cases]do
 if(ci-1)mod 3 ne SHARD then continue;end if;id:=cases[ci];q:=id eq 65 select 11 else 5;tag:=Sprintf("direct_case%o_A%o",id,q);md:=eval Read(root cat tag cat "_model_snapshot.m");b:=eval Read(root cat tag cat "_basis.m");res:=eval Read(root cat tag cat "_result.m");
 function DR(v)return B!(PK!v[1])/(B!(PK!v[2]));end function;
 f:=R![DR(v):v in md[1]];assert IsIrreducible(f);F<u>:=FunctionField(f);assert Genus(F)eq 25;OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
 function DF(v)return F![DR(a):a in v];end function;
 function Prime(pd)O:=pd[2]eq 1 select OF else OI;I:=ideal<O|[DF(v):v in pd[3]]>;assert IsPrime(I);P:=Place(I);assert Degree(P)eq pd[1];return P;end function;
 A0:=&+[DivisorGroup(F)|r[2]*Divisor(Prime(r[1])):r in b[6]];assert Degree(A0)eq b[1];BB:=DF(md[2][1]);half:=&+[Integers()!Floor(Rationals()!Valuation(BB,p)/2)*Divisor(p):p in Zeros(BB)cat Poles(BB)];assert Degree(half)eq -12;
 dims:=[];ranks:=[];
 for part in[0,1]do
  D:=A0+part*half;bf:=[DF(v):v in b[4+part]];n:=b[2+part];assert #bf eq n and &and[g ne 0:g in bf];I,J:=Ideals(D);
  for g in bf do assert IsIntegral(I*(g*OF))and IsIntegral(J*(g*OI));end for;
  den:=LCM(&cat[[Denominator(c):c in Eltseq(g)]:g in bf]);vv:=[[Numerator(c*(B!den)):c in Eltseq(g)]:g in bf];assert &and[&and[Denominator(c*(B!den))eq 1:c in Eltseq(g)]:g in bf];maxdeg:=Max(&cat[[Degree(c):c in v]:v in vv]);M:=Matrix(k,[&cat[[Coefficient(v[j],i):i in[0..maxdeg]]:j in[1..Degree(F)]]:v in vv]);rank:=Rank(M);assert rank eq n;
  dim:=Degree(D)-24;if Degree(D)le 48 then dim+:=#Basis(CanonicalDivisor(F)-D:Reduction:=false,Simplification:="None");end if;assert dim eq n;Append(~dims,dim);Append(~ranks,rank);
 end for;
 assert res[3]eq b[1]and res[4..5]eq dims and res[8]eq 0;Append(~out,[id,q,Degree(A0),dims[1],dims[2],ranks[1],ranks[2]]);print "C56_PARENT_BASE_BASIS_CASE",id,q,dims,ranks,Cputime(tt);
end for;Write(root cat Sprintf("parent_basis_poles_v2_s%o.m",SHARD),Sprint(out):Overwrite:=true);print "D31B_C56_PARENT_BASE_BASIS_COMPLETE",SHARD,#out,Cputime(tt);quit;
