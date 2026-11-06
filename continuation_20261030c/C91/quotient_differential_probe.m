load "work/deg10/continuation_20261030c/C91/quotient_model_common.m";

Omega,om:=SpaceOfHolomorphicDifferentials(CY); B:=[om(v):v in Basis(Omega)];
CM:=Matrix(GF(2),[Eltseq(Cartier(w)@@om):w in B]);
eta:=Differential(CY!K.1); hs:=[w/eta:w in B];
h:=hom<CY->CY|map<AF->CY|x:->CY!x>,zz+1>;
MH:=Matrix(GF(2),[Eltseq((h(q)*eta)@@om):q in hs]); I25:=IdentityMatrix(GF(2),25);
assert MH^2 eq I25 and MH*CM eq CM*MH;
U:=RowSpace((CM-I25)*(MH-I25));
printf "C91_DIFFERENTIALS_READY genus=%o CartierRank=%o hminusRank=%o Udim=%o time=%o\n",Genus(CY),Rank(CM),Rank(MH-I25),Dimension(U),Cputime(tt);
Write("work/deg10/continuation_20261030c/C91/differential_matrices.m",Sprintf("CM:=Matrix(GF(2),25,25,%o); MH:=Matrix(GF(2),25,25,%o); UB:=Matrix(GF(2),%o,25,%o);",Eltseq(CM),Eltseq(MH),Dimension(U),Eltseq(BasisMatrix(U))));
print "C91_DIFFERENTIAL_PROBE_COMPLETE"; quit;
