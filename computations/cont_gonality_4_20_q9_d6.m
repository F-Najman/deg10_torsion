SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(8*10^9);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
load "work/deg10/continuation_20261029/deps/dt_fast_differential_expansions.m";
F:=GF(9);X:=MDX11(4,20,F);C:=Curve(X);FF:=FunctionField(C);AF,af:=AlgorithmicFunctionField(FF);assert Genus(AF) eq 25;
function AP(x)
 pp,ee:=Support(FunctionFieldDivisor(Divisor(x)));assert #pp eq 1 and ee eq [1];return pp[1];
end function;
cs:=Cusps(X);rat:=[AP(c):c in cs|Degree(c) eq 1];quad:=[AP(c):c in cs|Degree(c) eq 2] cat [AP(c):c in NoncuspidalPlaces(X,2)];cubic:=[AP(c):c in NoncuspidalPlaces(X,3)];
assert #rat eq 24 and #quad eq 108 and #cubic eq 320;
O,om:=SpaceOfHolomorphicDifferentials(AF);ww:=[om(b):b in Basis(O)];assert #ww eq 25;
rjets:=[DifferentialExpansionMatrices(ww,p,4):p in rat];qjets:=[DifferentialExpansionMatrices(ww,p,1)[1]:p in quad];cjets:=[DifferentialExpansionMatrices(ww,p,1)[1]:p in cubic];
for i in [1..24] do for k in [1..4] do assert Dimension(RiemannRochSpace(k*Divisor(rat[i]))) eq k+1-Rank(HorizontalJoin(rjets[i][1..k]));end for;end for;
for i in [1..108] do assert Dimension(RiemannRochSpace(Divisor(quad[i]))) eq 3-Rank(qjets[i]);end for;
for i in [1..320] do assert Dimension(RiemannRochSpace(Divisor(cubic[i]))) eq 4-Rank(cjets[i]);end for;
function Enc(z)
 cc:=Eltseq(z,GF(3));assert #cc eq 2;return Integers()!cc[1]+3*(Integers()!cc[2]);
end function;
function EncMat(M)
 return [[Enc(a):a in Eltseq(row)]:row in Rows(M)];
end function;
records:=[* *];tested:=0;bad:=[* *];
procedure Check(A,B,mults,descr,~records,~tested,~bad)
 assert Degree(A) eq 3 and Nrows(B) eq 25 and Ncols(B) eq 3;
 U:=Nullspace(B);assert Dimension(U) eq 22;N:=Matrix(Basis(U));assert IsZero(N*B);
 J:=HorizontalJoin([rjets[i][mults[i]+1]:i in [1..24]]);M:=N*J;MT:=Transpose(M);
 for i in [1..22] do for j in [i+1..23] do for k in [j+1..24] do
  rank:=Rank(Submatrix(MT,[i,j,k],[1..22]));tested+:=1;
  if rank lt 3 then
   T:=A+Divisor(rat[i])+Divisor(rat[j])+Divisor(rat[k]);ell:=Dimension(RiemannRochSpace(T));assert ell eq 4-rank and ell ge 2;
   Append(~bad,[*descr,[i,j,k],ell*]);print "DECIC420_GONALITY6_FOUND",descr,[i,j,k],ell,tested,Cputime(tt);return;
  end if;
 end for;end for;end for;
 Append(~records,[*descr,EncMat(M)*]);
 if #records mod 100 eq 0 then print "DECIC420_GONALITY6_PROGRESS",#records,tested,Cputime(tt);end if;
end procedure;
for i in [1..24] do for j in [i..24] do for k in [j..24] do
 ms:=[0:l in [1..24]];ms[i]+:=1;ms[j]+:=1;ms[k]+:=1;A:=Divisor(rat[i])+Divisor(rat[j])+Divisor(rat[k]);B:=HorizontalJoin([HorizontalJoin(rjets[l][1..ms[l]]):l in [1..24]|ms[l] gt 0]);
 Check(A,B,ms,[*1,i,j,k*],~records,~tested,~bad);if #bad gt 0 then break;end if;
 end for;if #bad gt 0 then break;end if;end for;if #bad gt 0 then break;end if;end for;
if #bad eq 0 then for i in [1..108] do for j in [1..24] do
 ms:=[0:l in [1..24]];ms[j]:=1;A:=Divisor(quad[i])+Divisor(rat[j]);B:=HorizontalJoin(qjets[i],rjets[j][1]);
 Check(A,B,ms,[*2,i,j*],~records,~tested,~bad);if #bad gt 0 then break;end if;
end for;if #bad gt 0 then break;end if;end for;end if;
if #bad eq 0 then for i in [1..320] do
 Check(Divisor(cubic[i]),cjets[i],[0:l in [1..24]],[*3,i*],~records,~tested,~bad);if #bad gt 0 then break;end if;
end for;end if;
Write("work/deg10/continuation_20261029/M4N20_gonality6_arithmetic.txt",Sprint([*Coefficients(MinimalPolynomial(F.1)),records,tested,bad,[[EncMat(M):M in es]:es in rjets],[EncMat(M):M in qjets],[EncMat(M):M in cjets]*]));
print "DECIC420_GONALITY6_COMPLETE",#records,tested,bad,Cputime(tt);quit;
