SetColumns(0); SetQuitOnError(true); SetNthreads(1); SetMemoryLimit(12*10^9); tt:=Cputime();
F:=GF(2); K<t>:=FunctionField(F); P<v>:=PolynomialRing(K);
jt:=(t^2+5*t+13)*(t^4+7*t^3+20*t^2+19*t+1)^3/t;
f:=v^8+jt*v^7+v^4+1; assert IsIrreducible(f); AF<vv>:=FunctionField(f);
assert Genus(AF) eq 7; jaf:=AF!jt;
cs:=Poles(jaf); widths:=[-Valuation(jaf,c):c in cs];
assert #cs eq 4 and &and[Degree(c) eq 1:c in cs] and Seqset(widths) eq {1,7,13,91};
J:=JZero(91); fj:=FrobeniusPolynomial(J,2); if Type(fj) eq SeqEnum then fj:=&*fj; end if; nj:=Integers()!Evaluate(fj,1);
printf "C91_BASE_READY genus=%o widths=%o Jorder=%o time=%o\n",Genus(AF),widths,nj,Cputime(tt);
DG:=DivisorGroup(AF);Dzero:=DG!0;G,gm:=RayClassGroup(Dzero);
inv:=Invariants(G);assert &*[a:a in inv|a ne 0] eq nj and #[a:a in inv|a eq 0] eq 1;
printf "C91_PICARD invariants=%o time=%o\n",inv,Cputime(tt);
cols:=[];
for i in [1..4] do
 Append(~cols,Eltseq(Divisor(cs[i])@@gm));
 printf "C91_CUSP_CLASS index=%o class=%o time=%o\n",i,cols[i],Cputime(tt);
end for;
expect:=[w in {91,13} select 0 else 3:w in widths];
chars:=[];
for r in CartesianPower({0..3},Ngens(G)) do
 rr:=[r[i]:i in [1..Ngens(G)]];
 if not &and[inv[i]*rr[i] mod 4 eq 0:i in [1..#inv]] then continue;end if;
 if &and[(&+[cols[i][j]*rr[j]:j in [1..#inv]]) mod 4 eq expect[i]:i in [1..4]] then Append(~chars,rr);end if;
end for;
printf "C91_C4_CUSP_CHARACTERS count=%o images=%o time=%o\n",#chars,chars,Cputime(tt);

assert #chars eq 1;
C4:=AbelianGroup([4]);hh:=hom<G->C4|[C4![a]:a in chars[1]]>;
assert #Image(hh) eq 4;
fields,auts,wits:=InternalRayClassField(Dzero,Kernel(hh):WithAut:=true);
assert #fields eq #auts and #fields eq #wits and #fields eq 1;
wt:=wits[1][1];wc:=Eltseq(wt);WT:=wits[1][2];CY:=WT;E0:=BaseField(CY);
assert #wc eq 2 and Degree(CY) eq 2 and Degree(E0) eq 2 and BaseField(E0) cmpeq AF;
assert BaseRing(Parent(wt)) cmpeq AF;
Z0<z0>:=PolynomialRing(AF);Z1<z1>:=PolynomialRing(E0);
assert DefiningPolynomial(E0) eq z0^2+z0+wc[1];
assert DefiningPolynomial(CY) eq z1^2+z1+(E0!wc[2])+E0.1*(E0!wc[1]);
sigma0:=hom<E0->E0|map<AF->E0|x:->E0!x>,E0.1+1>;
sigma:=hom<CY->CY|map<E0->CY|x:->CY!sigma0(x)>,CY.1+CY!E0.1>;
assert sigma(sigma(sigma(sigma(CY.1)))) eq CY.1;
assert sigma(sigma(CY.1)) eq CY.1+1;
assert &and[sigma(CY!v) eq CY!v:v in [AF.1,AF!K.1,AF!K!BaseField(K).1]];
printf "C91_C4_CONSTRUCTED degree4=true cyclic4_verified=true time=%o\n",Cputime(tt);


assert Genus(CY) eq 25;
print "C91_WITT_TOWER_GENUS_VERIFIED";
Omega,om:=SpaceOfHolomorphicDifferentials(CY); B:=[om(v):v in Basis(Omega)];
CM:=Matrix(GF(2),[Eltseq(Cartier(w)@@om):w in B]);
eta:=Differential(CY!AF!K.1); hs:=[w/eta:w in B];
h:=hom<CY->CY|map<E0->CY|x:->CY!x>,CY.1+1>;
MH:=Matrix(GF(2),[Eltseq((h(q)*eta)@@om):q in hs]); I25:=IdentityMatrix(GF(2),25);
assert MH^2 eq I25 and MH*CM eq CM*MH;
U:=RowSpace((CM-I25)*(MH-I25));
printf "C91_WITT_DIFFERENTIALS_READY genus=%o CartierRank=%o hminusRank=%o Udim=%o time=%o\n",Genus(CY),Rank(CM),Rank(MH-I25),Dimension(U),Cputime(tt);
Write("work/deg10/continuation_20261030c/C91/witt_differential_matrices.m",Sprintf("CM:=Matrix(GF(2),25,25,%o); MH:=Matrix(GF(2),25,25,%o); UB:=Matrix(GF(2),%o,25,%o);",Eltseq(CM),Eltseq(MH),Dimension(U),Eltseq(BasisMatrix(U))));
print "C91_WITT_DIFFERENTIAL_PROBE_COMPLETE"; quit;
