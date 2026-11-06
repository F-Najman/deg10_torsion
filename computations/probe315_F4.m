SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(3*10^9);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
X:=MDX11(3,15,GF(4));C:=Curve(X);FF:=FunctionField(C);AF,afmap:=AlgorithmicFunctionField(FF);g:=Genus(AF);assert g eq 9;
O,om:=SpaceOfHolomorphicDifferentials(AF);B:=[om(v):v in Basis(O)];F:=GF(4);
CM:=Matrix(F,[Eltseq(Cartier(Cartier(w)) @@ om):w in B]);U:=RowSpace(CM-IdentityMatrix(F,g));UB:=BasisMatrix(U);print "DECIC_315_F4_CARTIER_SQUARE_IMAGE",Dimension(U),Cputime(tt);
jaf:=afmap(jInvariant(EllipticCurve(X`_E)));cs:=Poles(jaf);ds:=[Degree(c):c in cs];assert {*d:d in ds*} eq {*1^^16,2^^8*};
VV:=VectorSpace(F,g);rows:=[* *];
code:=func<a|Integers()!Eltseq(a,GF(2))[1]+2*(Integers()!Eltseq(a,GF(2))[2])>;
for i in [1..#cs] do
 for m in [1..5 div ds[i]] do
  V,vm:=DifferentialSpace(m*Divisor(cs[i]));W:=sub<VV|[Eltseq(vm(v) @@ om):v in Basis(V)]>;
  assert Dimension(W) eq g-m*ds[i];J:=Transpose(BasisMatrix(Nullspace(Transpose(BasisMatrix(W)))));assert Nrows(J) eq g and Ncols(J) eq m*ds[i];A:=UB*J;
  Append(~rows,[*i,m,Nrows(A),Ncols(A),[code(a):a in Eltseq(A)]*]);
 end for;
 print "DECIC_315_F4_RESTRICTIONS",i,#cs,Cputime(tt);
end for;
Write("work/deg10/notes/probe315_F4_data.txt",Sprint([*ds,Dimension(U),rows*]));
print "DECIC_315_F4_DIAGNOSTIC_EXPORT_COMPLETE",Cputime(tt);quit;
