load "work/deg10/continuation_20261030b/C77/geometry_common.m";
SetNthreads(1);t0:=Cputime();Om,om:=SpaceOfHolomorphicDifferentials(CY);bs:=[om(v):v in Basis(Om)];C:=Matrix(GF(2),[Eltseq(Cartier(w)@@om):w in bs]);U:=RowSpace(C*(C-IdentityMatrix(GF(2),37)));assert Dimension(U) eq 32;
rat:=[cys[i]:i in [1..#cys]|cywidths[i] eq 77 and Degree(cys[i]) eq 1][1];quad:=[cys[i]:i in [1..#cys]|cywidths[i] eq 7 and Degree(cys[i]) eq 2][1];out:=[];
for D in [10*Divisor(rat),5*Divisor(quad)] do
 V,vm:=DifferentialSpace(D);W:=sub<VectorSpace(GF(2),37)|[Eltseq(vm(w)@@om):w in Basis(V)]>;rank:=Dimension(U)-Dimension(U meet W);assert rank eq 10;Append(~out,[Degree(D),Dimension(W),rank]);print"C77_FI_NATIVE_CONTROL",#out,rank,Cputime(t0);
end for;
Write("work/deg10/continuation_20261030b/C77/fi_native_controls.m",Sprint(out):Overwrite:=true);print"C77_FI_NATIVE_CONTROLS_COMPLETE",#out,Cputime(t0);quit;
