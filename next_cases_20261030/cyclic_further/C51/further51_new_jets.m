SetColumns(0);SetQuitOnError(true);SetMemoryLimit(16*10^9);tt:=Cputime();
AttachSpec("work/degree9/shared/mdmagma/v2/mdmagma.spec");
X:=MDX1(51,GF(2));AF,am:=AlgorithmicFunctionField(FunctionField(Curve(X)));assert Genus(AF) eq 65;
Omega,om:=SpaceOfHolomorphicDifferentials(AF);B:=[om(v):v in Basis(Omega)];CMnew:=Matrix(GF(2),[Eltseq(Cartier(w)@@om):w in B]);Unew:=RowSpace(CMnew-IdentityMatrix(GF(2),65));assert Dimension(Unew) eq 62;
j:=am(jInvariant(EllipticCurve(X`_E)));cs:=Poles(j);assert Degree(j) eq 1152 and {*Degree(c):c in cs*} eq {*1^^16,2^^8,8^^4*};rats:=[c:c in cs|Degree(c) eq 1];quads:=[c:c in cs|Degree(c) eq 2];
printf "DECIC51_NEW_JET_SETUP genus=65 Udim=62 rational=16 quadratic=8 octic=4 time=%o\n",Cputime(tt);
eta:=Differential(SeparatingElement(AF));ratios:=[w/eta:w in B];
function CompletionRestriction(AF,ratios,eta,P,m)
 k:=Valuation(eta,P);prec:=m+Maximum(0,-k)+2;
 repeat
  L,ex:=Completion(AF,P:Precision:=prec);t:=L.1;es:=[t^k*ex(f):f in ratios];good:=&and[AbsolutePrecision(s) ge m:s in es];
  if not good then prec*:=2;end if;
 until good;
 assert &and[IsZero(s) or Valuation(s) ge 0:s in es];
 A:=Matrix(GF(2),[&cat[Eltseq(Coefficient(s,j),GF(2)):j in [0..m-1]]:s in es]);
 assert Nrows(A) eq 65 and Ncols(A) eq m*Degree(P);return A;
end function;
AR10:=CompletionRestriction(AF,ratios,eta,rats[1],10);assert Rank(BasisMatrix(Unew)*AR10) eq 10;
printf "DECIC51_NEW_RATIONAL_CAP10_DONE rank10=10 time=%o\n",Cputime(tt);
AQ5:=CompletionRestriction(AF,ratios,eta,quads[1],5);assert Rank(BasisMatrix(Unew)*AQ5) eq 10;
printf "DECIC51_NEW_QUADRATIC_CAP5_DONE rank10=10 time=%o\n",Cputime(tt);
load "work/degree9/cyclic/continuation_20261026/N51_geometry/cyclic_sieves/N51_frozen_jets.m";
n:=65;VO:=VectorSpace(GF(2),n*n);eqs:=[];Wold8:=Nullspace(basisjets[1]);Anew8:=ColumnSubmatrix(AR10,8);Wnew8:=Nullspace(Anew8);assert Dimension(Wold8) eq 57 and Dimension(Wnew8) eq 57;
for i,j in [1..n] do
 row:=VO!0;for k in [1..n] do row[(k-1)*n+j]+:=CartierMatrix[i,k];row[(i-1)*n+k]+:=CMnew[k,j];end for;Append(~eqs,row);
end for;
for w in Basis(Wold8) do for c in [1..8] do
 row:=VO!0;for i,j in [1..n] do row[(i-1)*n+j]:=w[i]*Anew8[j,c];end for;Append(~eqs,row);
end for;end for;
Ktrans:=Nullspace(Transpose(Matrix(eqs)));assert Dimension(Ktrans) eq 1;Ttrans:=Matrix(GF(2),n,n,Eltseq(Basis(Ktrans)[1]));assert Rank(Ttrans) eq 65;
assert CartierMatrix*Ttrans eq Ttrans*CMnew and RowSpace(BasisMatrix(Wold8)*Ttrans) eq Wnew8 and RowSpace(ImageBasis*Ttrans) eq Unew;
printf "DECIC51_BASIS_ALIGNMENT_DONE dimension=1 invertible=true anchor_degree=8 time=%o\n",Cputime(tt);
ratindices:=[i:i in [1..#ds]|ds[i] eq 1];assert #ratindices eq 16 and ratindices[1] eq 1;
ratA:=[Ttrans^-1*ColumnSubmatrix(basisjets[i],1):i in ratindices];R9:=ColumnSubmatrix(AR10,9);checked:=0;
for k in [2..16] do assert Rank(BasisMatrix(Unew)*HorizontalJoin(R9,ratA[k])) eq 10;checked+:=1;end for;
assert checked eq 15;
printf "DECIC51_NINE_PLUS_ONE_RANKS_DONE cases=15 rank=10 time=%o\n",Cputime(tt);
f:="work/deg10/next_cases_20261030/cyclic_further/C51/new_jets.json";out:=Open(f,"w");
fprintf out,"{\"level\":51,\"degree\":10,\"genus\":65,\"U_dimension\":62,\"Cartier\":%o,\"U_basis\":%o,\"T_old_to_new\":%o,\"anchor_kernel_new\":%o,\"rational_cap10\":%o,\"quadratic_cap5\":%o,\"rational_first_jets\":%o,\"rational_old_indices\":%o,\"alignment_dimension\":1,\"rational_cap10_rank\":10,\"quadratic_cap5_rank\":10,\"nine_plus_one_cases\":15,\"complete\":true}\n",[Integers()!x:x in Eltseq(CMnew)],[Integers()!x:x in Eltseq(BasisMatrix(Unew))],[Integers()!x:x in Eltseq(Ttrans)],[Integers()!x:x in Eltseq(BasisMatrix(Wnew8))],[Integers()!x:x in Eltseq(AR10)],[Integers()!x:x in Eltseq(AQ5)],[[Integers()!x:x in Eltseq(a)]:a in ratA],ratindices;delete out;
printf "DECIC51_NEW_JETS_COMPLETE rational_cap10=true quadratic_cap5=true nine_plus_one=15 unique_basis_alignment=true time=%o\n",Cputime(tt);quit;
