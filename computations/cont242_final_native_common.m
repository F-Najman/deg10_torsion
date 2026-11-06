SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(6*10^9);tt:=Cputime();
root:="work/deg10/next_cases_20261030/C242/direct_hecke17/";
a:=eval Read(root cat TAG cat "/model_bases.m");k:=GF(5);PP:=PolynomialRing(k);K<x>:=FunctionField(k);R<Y>:=PolynomialRing(K);
function DR(v) return K!(PP!v[1])/K!(PP!v[2]);end function;
F<y>:=FunctionField(R![DR(v):v in a[2]]);assert Degree(F) eq a[1] and Genus(F) eq 25;
function DF(v) return F![DR(t):t in v];end function;
OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
function Mat(bs) return Matrix(K,[Eltseq(F!v):v in bs]);end function;
function IntegralEntries(M,inf)
 return &and[inf select Degree(Numerator(v)) le Degree(Denominator(v)) else Degree(Denominator(v)) eq 0:v in Eltseq(M)];
end function;
procedure CheckLattice(bs,ns,inf)
 assert #bs eq Degree(F) and #ns eq Degree(F);
 B:=Mat(bs);N:=Mat(ns);assert Determinant(B) ne 0 and Determinant(N) ne 0;
 T:=B*N^-1;assert IntegralEntries(T,inf) and IntegralEntries(T^-1,inf);
end procedure;
BF:=[DF(v):v in a[4]];BI:=[DF(v):v in a[5]];
CheckLattice(BF,[F!v:v in Basis(OF)],false);CheckLattice(BI,[F!v:v in Basis(OI)],true);
rows:=eval Read(root cat "selected_native_replay_cases.m");out:=[* *];
for row in rows do
 if row[1] ne TAG then continue;end if;
 ci:=row[2];mi:=row[3];ex:=row[4];assert ex eq 0;
 b:=eval Read(root cat TAG cat Sprintf("/ideal_basis_case%o_map%o.m",ci,mi));assert b[1] eq ci and b[2] eq mi;
 ib:=[DF(v):v in b[3]];jb:=[DF(v):v in b[4]];I:=ideal<OF|ib>;J:=ideal<OI|jb>;
 CheckLattice(ib,[F!v:v in Basis(I)],false);CheckLattice(jb,[F!v:v in Basis(J)],true);
 D:=-Divisor(I,J);assert Degree(D) eq 0;ID,JD:=Ideals(D);assert ID eq I^-1 and JD eq J^-1;
 hh:=#Basis(D:Reduction:=false,Simplification:="None");assert hh eq ex;
 Append(~out,[*ci,mi,hh*]);print "C242_FINAL_SECOND_NATIVE_TEST",TAG,ci,mi,hh,Cputime(tt);
end for;
assert #out eq (TAG eq "A11pure" select 5 else 2);
Write(root cat TAG cat "/parent_native_RR_replay.m",Sprint(out):Overwrite:=true);
print "C242_FINAL_SECOND_NATIVE_COMPLETE",TAG,#out,Cputime(tt);quit;
