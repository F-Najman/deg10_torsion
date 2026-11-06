SetColumns(0);SetQuitOnError(true);SetMemoryLimit(12*10^9);
AttachSpec("work/degree9/shared/mdmagma/v2/mdmagma.spec");
load "work/deg10/next_cases_20261030/cyclic_primepowers/deps/dt_expansions_local.m";
tt:=Cputime();X:=MDX1(49,GF(2));AF,am:=AlgorithmicFunctionField(FunctionField(Curve(X)));assert Genus(AF) eq 69;
b:=am(X`_coordinates[1]);c:=am(X`_coordinates[2]);j:=(1+c)^12/(b^3*(b+c*(c+1)^3));
Om,om:=SpaceOfHolomorphicDifferentials(AF);B:=[om(v):v in Basis(Om)];assert #B eq 69;
cs:=Poles(j);assert Degree(j) eq 1176 and {*Degree(P):P in cs*} eq {*1^^21,3^^6,21*};rat:=[P:P in cs|Degree(P) eq 1];cub:=[P:P in cs|Degree(P) eq 3];
printf "DECIC_C49_COMPLETION_CONTROL_SETUP genus=69 time=%o\n",Cputime(tt);
function NativeCompletionMatrix(B,P,m)
 AF:=FunctionField(P);du:=UniformizingDifferential(B,P);L,lc:=Completion(AF,P:Precision:=m+2);rows:=[];
 for w in B do
  ser:=lc(AF!RationalFunction(w/du));assert Valuation(ser) ge 0;
  Append(~rows,&cat[Eltseq(Coefficient(ser,k),GF(2)):k in [0..m-1]]);
 end for;
 return Matrix(GF(2),rows);
end function;
A:=HorizontalJoin(DifferentialExpansionMatrices(B,rat[1],7));NA:=NativeCompletionMatrix(B,rat[1],7);
assert Nrows(A) eq 69 and Ncols(A) eq 7 and Nrows(NA) eq 69 and Ncols(NA) eq 7;
for m in [1..7] do assert Nullspace(Submatrix(A,1,1,69,m)) eq Nullspace(Submatrix(NA,1,1,69,m));end for;
printf "DECIC_C49_COMPLETION_RATIONAL_ALL_SEVEN_KERNELS_AGREE time=%o\n",Cputime(tt);
Q:=HorizontalJoin(DifferentialExpansionMatrices(B,cub[1],1));NQ:=NativeCompletionMatrix(B,cub[1],1);
assert Nrows(Q) eq 69 and Ncols(Q) eq 3 and Nullspace(Q) eq Nullspace(NQ);
AA:=HorizontalJoin(A,Q);NN:=HorizontalJoin(NA,NQ);assert Nullspace(AA) eq Nullspace(NN) and Rank(AA) eq 10;
printf "DECIC_C49_COMPLETION_CONTROL_COMPLETE degree=10 full_kernel_agrees=true kernel_dim=%o direct_rank=10 built_in_completion=true time=%o\n",Dimension(Nullspace(NN)),Cputime(tt);quit;
