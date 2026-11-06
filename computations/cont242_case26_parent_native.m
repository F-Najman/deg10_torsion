SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(4*10^9);tt:=Cputime();
root:="work/deg10/next_cases_20261030/C242/geometry/";m:=eval Read(root cat "case26_section_model.m");a:=eval Read(root cat "case26_section_RR.m");assert #a[2] eq 0;
k5:=GF(5);P5:=PolynomialRing(k5);k<zeta>:=ext<k5|P5!m[1]>;PP:=PolynomialRing(k);K<x>:=FunctionField(k);R<Y>:=PolynomialRing(K);
function EC(v)return Evaluate(P5!v,zeta);end function;
function DR(v)return K!(PP![EC(t):t in v[1]])/K!(PP![EC(t):t in v[2]]);end function;
F<y>:=FunctionField(R![DR(v):v in m[2]]);assert Genus(F) eq 17 and Degree(F) eq 8;
function DF(v)return F![DR(t):t in v];end function;
OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
function ED(rows)
 D:=DivisorGroup(F)!0;
 for row in rows do p,e:=Explode(row);O:=p[1] eq 1 select OF else OI;I:=ideal<O|[DF(t):t in p[3]]>;assert IsPrime(I);Q:=Place(I);assert Degree(Q) eq p[2];D+:=e*Divisor(Q);end for;return D;
end function;
function Mat(bs)return Matrix(K,[Eltseq(F!b):b in bs]);end function;
function IntEntries(T,inf)return &and[inf select Degree(Numerator(t)) le Degree(Denominator(t)) else Degree(Denominator(t)) eq 0:t in Eltseq(T)];end function;
procedure Lattice(bs,ns,inf)
 assert #bs eq 8 and #ns eq 8;B:=Mat(bs);N:=Mat(ns);assert Determinant(B) ne 0 and Determinant(N) ne 0;T:=B*N^-1;assert IntEntries(T,inf) and IntEntries(T^-1,inf);
end procedure;
u:=DF(m[4]);fk:=DF(m[5]);s:=u^2+8;D0:=s^2-3*s+9;A0:=2*s^3-9*s^2+27*s-27;assert A0^2-4*D0^3 eq -27*D0^2;zs:=Roots(PolynomialRing(k)![1,1,1]);assert #zs eq 2 and &or[fk eq (A0+3*(t[1]-t[1]^2)*D0)/2:t in zs];
ps,es:=Support(Divisor(fk));assert &+[Degree(ps[i]):i in [1..#ps]|es[i] mod 3 ne 0] eq 12;
chars:=[ED(v):v in m[6]];assert chars[1] eq DivisorGroup(F)!0;
for j in [1,2] do want:=&+[Floor(j*es[i]/3)*Divisor(ps[i]):i in [1..#ps]];assert chars[j+1] eq want;end for;
assert [Degree(D):D in chars] eq [0,-6,-6];rats:=[ED(v):v in m[7]];cubs:=[ED(v):v in m[8]];assert [Degree(D):D in rats] eq [1:i in [1..6]] and [Degree(D):D in cubs] eq [3,3];qs:=AssociativeArray();for r in m[9] do qs[r[1]]:=ED(r[2]);assert Degree(qs[r[1]]) eq 2;end for;
out:=[* *];n:=0;
for row in a[1] do
 qi,ci,ri,dims,proofs:=Explode(row);AA:=qs[qi]+cubs[ci]+rats[ri];assert Degree(AA) eq 6 and dims eq [1,0,0];
 for j in [1..3] do
  p:=proofs[j];D:=ED(p[2]);assert D eq AA+chars[j];I,J:=Ideals(D);ib:=[DF(v):v in p[3]];jb:=[DF(v):v in p[4]];assert ideal<OF|ib> eq I and ideal<OI|jb> eq J;Lattice(ib,[F!v:v in Basis(I)],false);Lattice(jb,[F!v:v in Basis(J)],true);
  bb:=Basis(D:Reduction:=false,Simplification:="None");assert #bb eq p[1] and #bb eq dims[j];assert #p[5] eq #bb;
  for v in p[5] do h:=DF(v);assert h ne 0 and IsEffective(D+Divisor(h));end for;
  n+:=1;Append(~out,[*qi,ci,ri,j,#bb*]);
 end for;
 print "C242_CASE26_PARENT_NATIVE_PATTERN",qi,ci,ri,n,Cputime(tt);
end for;assert #a[1] eq 24 and n eq 72;
Write(root cat "case26_parent_native_replay.m",Sprint(out):Overwrite:=true);print "C242_CASE26_PARENT_NATIVE_COMPLETE",24,n,Cputime(tt);quit;
