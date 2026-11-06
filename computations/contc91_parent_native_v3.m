SetColumns(0);SetNthreads(1);SetQuitOnError(true);SetMemoryLimit(8*10^9);Alarm(500);tt:=Cputime();
root:="work/deg10/continuation_20261030c/C91/";load "work/deg10/continuation_20261030c/C91/native_packet.m";PCM:=CM;PMT:=MT;PU:=RowSpace(UB);
k:=GF(2);PP:=PolynomialRing(k);K<t>:=FunctionField(k);R<X>:=PolynomialRing(K);
function DR(a)return K!(PP!a[1])/K!(PP!a[2]);end function;
f:=R![DR(v):v in DefPolyER];assert Degree(f) eq 14 and IsIrreducible(f);Y<s>:=FunctionField(f);assert Genus(Y) eq 23;
function DF(a)return Y![DR(v):v in a];end function;
b:=t^3+t^2;c:=t^2+t;E:=EllipticCurve([K|1-c,b,b,0,0]);P:=E![0,0];assert Order(P) eq 7;j7:=Y!jInvariant(E);
j13:=(s^2+s+1)*(s^4+s^3+s+1)^3/s;assert j13 eq j7;
z:=1/s;j13dual:=(z^2+z+1)*(z^4+z^3+z+1)^3/z;Phi:=ChangeRing(ClassicalModularPolynomial(13),k);assert Evaluate(Phi,[j13,j13dual]) eq 0;
tau:=hom<Y->Y|hom<K->Y|1/(Y!t+1)>,s>;assert tau(tau(tau(Y!t))) eq Y!t and tau(tau(tau(s))) eq s and tau(j7) eq j7;
O:=MaximalOrderFinite(Y);OI:=MaximalOrderInfinite(Y);primes:=[];
for i in[1..#PrimeER] do
 row:=PrimeER[i];ff:=DF(row[2]);gg:=DF(row[3]);cand:=[];
 for OO in[O,OI] do
  try
   aa:=OO!ff;bb:=OO!gg;I:=ideal<OO|aa,bb>;
   if IsPrime(I) then Q:=Place(I);if Degree(Q) eq row[1] then Append(~cand,Q);end if;end if;
  catch err
   skipped:=true;
  end try;
 end for;
 assert #cand eq 1;Q:=cand[1];assert Valuation(ff,Q) gt 0 and Valuation(gg,Q) gt 0;Append(~primes,Q);
 if i le 4 then assert -Valuation(j7,Q) in[13,91];else
  assert Degree(Q) eq 2 and Valuation(Y!(t^2+t+1),Q) gt 0 and Valuation(j7,Q) gt 0;
  u:=UniformizingElement(Q);assert Valuation(tau(u),Q) eq 1;lam:=Evaluate(tau(u)/u,Q);assert lam^3 eq 1 and lam ne 1;
 end if;
end for;assert #primes eq 6;
V,om:=SpaceOfHolomorphicDifferentials(Y);B:=[om(v):v in Basis(V)];eta:=Differential(Y!t);
CM:=Matrix(k,[Eltseq(Cartier(w)@@om):w in B]);MT:=Matrix(k,[Eltseq((tau(w/eta)*Differential(tau(Y!t)))@@om):w in B]);I:=IdentityMatrix(k,23);
assert CM eq PCM and MT eq PMT and MT^3 eq I and CM*MT eq MT*CM;U:=RowSpace((CM-I)*(MT-I));assert Dimension(U) eq 14 and U eq PU;
VO:=VectorSpace(k,23);kslocal:=AssociativeArray();n1:=0;n2:=0;cmranks:=[];
for rec in KernelRecords do
 kind,pi,m,ra,rb,dm,flat:=Explode(rec);D:=m*Divisor(primes[pi]);
 if kind eq 2 then D+:=Divisor(primes[ra])+Divisor(primes[rb]);assert Degree(D) eq 10;end if;
 Z,zm:=DifferentialSpace(D);W:=sub<VO|[VO!Eltseq(zm(v)@@om):v in Basis(Z)]>;
 saved:=RowSpace(Matrix(k,dm,23,flat));assert W eq saved and Dimension(W) eq 23-Degree(D);
 if kind eq 1 then kslocal[<pi,m>]:=W;n1+:=1;else
  rr:=Dimension(U)-Dimension(U meet W);assert rr eq 8;Append(~cmranks,[pi,ra,rb,rr]);n2+:=1;
 end if;
end for;assert n1 eq 33 and n2 eq 12;
rat:=[i:i in[1..4]|Degree(primes[i]) eq 1];cub:=[i:i in[1..4]|Degree(primes[i]) eq 3];assert #rat eq 3 and #cub eq 1;
ranks:=[];
for a in[0..10] do for b in[0..10-a] do for c in[0..10-a-b] do
 left:=10-a-b-c;if left mod 3 ne 0 then continue;end if;d:=left div 3;mult:=[0:i in[1..4]];mult[rat[1]]:=a;mult[rat[2]]:=b;mult[rat[3]]:=c;mult[cub[1]]:=d;
 W:=VO;for i in[1..4] do if mult[i] gt 0 then W:=W meet kslocal[<i,mult[i]>];end if;end for;
 kk:=Min([a,b,c])+d;rr:=Dimension(U)-Dimension(U meet W);assert rr eq 10-kk;Append(~ranks,[a,b,c,d,kk,rr]);
end for;end for;end for;assert #ranks eq 120;
Write(root cat "parent_native_replay.m",Sprint([*ranks,cmranks*]):Overwrite:=true);print "C91_PARENT_NATIVE_COMPLETE",#primes,n1+n2,#ranks,#cmranks,Cputime(tt);quit;
