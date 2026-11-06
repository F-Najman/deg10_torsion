SetColumns(0);SetQuitOnError(true);SetMemoryLimit(32*10^9);
AttachSpec("work/degree9/shared/mdmagma/v2/mdmagma.spec");load "work/deg10/next_cases_20261030/cyclic_primepowers/deps/dt_expansions_local.m";
tt:=Cputime();X:=MDX1(49,GF(2));FF:=FunctionField(Curve(X));AF,am:=AlgorithmicFunctionField(FF);assert Genus(AF) eq 69;
u:=am(X`_coordinates_xy[1]);v:=am(X`_coordinates_xy[2]);K:=BaseRing(AF);assert u eq AF!K.1 and v eq AF.1;
E:=EllipticCurve([am(z):z in X`_E]);P:=E![AF!0,AF!0];xy:=TateNormalForm_xy(E,2*P);sig0:=hom<K->AF|xy[1]>;dg:=hom<AF->AF|sig0,xy[2]>;
Om,om:=SpaceOfHolomorphicDifferentials(AF);B:=[om(z):z in Basis(Om)];CM:=Matrix(GF(2),[Eltseq(Cartier(w)@@om):w in B]);eta:=Differential(u);deta:=Differential(dg(u));MD:=Matrix(GF(2),[Eltseq((dg(w/eta)*deta)@@om):w in B]);I:=IdentityMatrix(GF(2),69);assert MD^21 eq I and MD^7 ne I and MD^3 ne I and MD*CM eq CM*MD;
U:=RowSpace(CM-I);assert Dimension(U) eq 68 and RowSpace(BasisMatrix(U)*MD) eq U;
printf "DECIC_C49_CARTIER_TRANSPORT_SETUP genus=69 Udim=68 diamond2_exact_order=21 time=%o\n",Cputime(tt);
b:=am(X`_coordinates[1]);c:=am(X`_coordinates[2]);j:=(1+c)^12/(b^3*(b+c*(c+1)^3));cs:=Poles(j);assert Degree(j) eq 1176 and {*Degree(z):z in cs*} eq {*1^^21,3^^6,21*};rats:=[z:z in cs|Degree(z) eq 1];cub:=[z:z in cs|Degree(z) eq 3];
f23:=[am(z):z in X`_coordinates_F2F3];sigs:=[[Valuation(z,pl):z in f23]:pl in cub];assert #Seqset(sigs) eq 3 and &and[Multiplicity(sigs,z) eq 2:z in Seqset(sigs)];pair:=[pl:pl in cub|[Valuation(z,pl):z in f23] eq sigs[1]];assert #pair eq 2;
// delta2 has two cubic closed-cusp orbits of length3; this is also checked
// on primitive Tate labels in the separate read-only arithmetic replay.
for pl in pair do
 ss:=[[Valuation(z,pl):z in f23],[Valuation(dg(z),pl):z in f23],[Valuation(dg(dg(z)),pl):z in f23]];
 assert #Seqset(ss) eq 3 and Seqset(ss) eq Seqset(sigs);
end for;
AA:=HorizontalJoin(DifferentialExpansionMatrices(B,rats[1],7));assert Nrows(AA) eq 69 and Ncols(AA) eq 7;
BB:=[HorizontalJoin(DifferentialExpansionMatrices(B,pl,3)):pl in pair];assert &and[Nrows(z) eq 69 and Ncols(z) eq 9:z in BB];
for z in BB do for m in [1..3] do W:=Nullspace(Submatrix(z,1,1,69,3*m));assert sub<Generic(W)|[w*MD^3:w in Basis(W)]> eq W;end for;end for;
ds:=[1:i in [1..19]]cat[3:i in [1..6]];caps:=[7:i in [1..19]]cat[3:i in [1..6]];
f:="work/deg10/next_cases_20261030/cyclic_primepowers/c49_cartier_jets.m";
PrintFile(f,Sprintf("N:=49;D:=10;gx:=69;udim:=68;ds:=%o;caps:=%o;basisjets:=AssociativeArray();localjets:=AssociativeArray();",ds,caps):Overwrite:=true);
PrintFile(f,Sprintf("CartierMatrix:=Matrix(GF(2),69,69,%o);DiamondMatrix:=Matrix(GF(2),69,69,%o);ImageBasis:=Matrix(GF(2),68,69,%o);",Eltseq(CM),Eltseq(MD),Eltseq(BasisMatrix(U))));
blocks:=[* MD^k*AA:k in [0..18] *];
for k in [0..2] do for z in BB do Append(~blocks,MD^k*z);end for;end for;
assert #blocks eq 25;
for i in [1..25] do
 z:=blocks[i];PrintFile(f,Sprintf("basisjets[%o]:=Matrix(GF(2),69,%o,%o);",i,Ncols(z),Eltseq(z)));
 for m in [1..caps[i]] do JJ:=BasisMatrix(U)*Submatrix(z,1,1,69,ds[i]*m);assert Rank(JJ) eq ds[i]*m;PrintFile(f,Sprintf("localjets[<%o,%o>]:=Matrix(GF(2),68,%o,%o);",i,m,Ncols(JJ),Eltseq(JJ)));end for;
end for;
PrintFile(f,"print \"DECIC_C49_TRANSPORT_MIXED10_JET_DATA_READY\";");
printf "DECIC_C49_CARTIER_TRANSPORT_EXPORT_COMPLETE common_basis=true selected_rational19=true cubic_orbits3plus3=true time=%o\n",Cputime(tt);quit;
