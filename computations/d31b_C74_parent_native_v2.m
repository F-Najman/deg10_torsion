SetColumns(0);SetNthreads(1);SetQuitOnError(true);SetMemoryLimit(10*10^9);Alarm(650);tt:=Cputime();AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");root:="work/deg10/continuation_20261031b/C74/";md:=eval Read(root cat"formal_model.m");cj:=eval Read(root cat"rational_cusp_jets.m");sx:=eval Read(root cat"sextic_restriction.m");decs:=eval Read(root cat"decic_restrictions.m");assert #decs eq 99;k:=GF(2);PK:=PolynomialRing(k);K<t>:=FunctionField(k);R<Z>:=PolynomialRing(K);
function DR(v)return K!(PK!v[1])/(K!(PK!v[2]));end function;
f:=R![DR(v):v in md[1]];assert IsIrreducible(f);F<u>:=FunctionField(f);assert Genus(F)eq 40;OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);eta:=Differential(F!t);
function DF(v)return F![DR(w):w in v];end function;
function DM(v)return Matrix(k,v[1],v[2],v[3]);end function;
function Prime(pd)O:=pd[2]eq 1 select OF else OI;gs:=[DF(v):v in pd[3]];assert #gs eq 2;I:=ideal<O|gs>;assert IsIntegral(I)and IsPrime(I);P:=Place(I);assert Degree(P)eq pd[1]and Ideal(P)eq I;return P;end function;
B:=[DF(v)*eta:v in md[2]];Om,om:=SpaceOfHolomorphicDifferentials(F);assert Dimension(Om)eq 40;T:=Matrix(k,[Eltseq(w@@om):w in B]);assert Rank(T)eq 40;TI:=T^-1;MD:=DM(md[3]);CM:=DM(md[4]);U:=RowSpace(DM(md[5]));assert Dimension(U)eq 38 and U eq RowSpace(MD-IdentityMatrix(k,40))and U eq RowSpace((CM-IdentityMatrix(k,40))*(MD-IdentityMatrix(k,40)));
xx:=F!t;yy:=F.1;r:=(xx^2*yy-xx*yy+yy-1)/(xx^2*yy-xx);s:=(xx*yy-yy+1)/(xx*yy);b:=r*s*(r-1);c:=s*(r-1);E:=EllipticCurve([F|1-c,-b,-b,0,0]);P0:=E![0,0];assert 37*P0 eq E!0 and P0 ne E!0;vv:=TateNormalForm_xy(E,2*P0);assert vv eq [DF(z):z in md[6]];h0:=hom<K->F|vv[1]>;h:=hom<F->F|h0,vv[2]>;de:=Differential(h(xx));
for i in[1..40]do assert h(B[i]/eta)*de eq &+[MD[i,j]*B[j]:j in[1..40]];assert Cartier(B[i])eq &+[CM[i,j]*B[j]:j in[1..40]];end for;
assert MD^18 eq 1 and MD*CM eq CM*MD;print "C74_PARENT_ACTUAL_DIAMOND_CARTIER_BASIS_PASS",Cputime(tt);
function Key2(vals)
 L:=Parent(vals[1]);d:=LCM([Degree(MinimalPolynomial(v,k)):v in vals]);ok:=false;
 for bits in[[1,0],[0,1],[1,1]]do g:=&+[L!bits[i]*vals[i]:i in[1,2]];mp:=MinimalPolynomial(g,k);if Degree(mp)eq d then chosen:=bits;ok:=true;break;end if;end for;assert ok;
 MM:=Matrix(k,[Eltseq(g^i,k):i in[0..d-1]]);assert Rank(MM)eq d;co:=[Eltseq(Solution(MM,Vector(k,Eltseq(v,k)))):v in vals];return[*d,chosen,[Integers()!v:v in Coefficients(mp)],[[Integers()!v:v in row]:row in co]*];
end function;
A0:=DM(cj[1]);rats:=[Prime(p):p in cj[4]];assert #rats eq 18 and #Seqset(rats)eq 18;P6:=Prime(sx[1]);A6:=DM(sx[2]);assert Degree(P6)eq 6;
function EvalRestriction(P)
 ii:=0;for j in[1..#B]do if Valuation(B[j],P)eq 0 then ii:=j;break;end if;end for;assert ii gt 0;ratios:=[F!RationalFunction(w/B[ii]):w in B];assert &and[Valuation(v,P)ge 0:v in ratios];AA:=Matrix(k,[Eltseq(Evaluate(v,P),k):v in ratios]);assert Nrows(AA)eq 40 and Ncols(AA)eq Degree(P);return AA,ratios;
end function;
Afirst,ratios:=EvalRestriction(rats[1]);L,ev:=Completion(F,rats[1]:Precision:=160);series:=[ev(v):v in ratios];assert &and[Valuation(v)ge 0 and AbsolutePrecision(v)ge 10:v in series];AC:=Matrix(k,[[k!Coefficient(v,j):j in[0..9]]:v in series]);
for m in[1..10]do assert Nullspace(ColumnSubmatrix(AC,1,m))eq Nullspace(ColumnSubmatrix(A0,1,m));end for;
for j in[1..18]do AA,rr:=EvalRestriction(rats[j]);power:=cj[2][j]-1;assert Nullspace(AA)eq Nullspace((MD^power)*ColumnSubmatrix(A0,1,1));end for;
AS,rr:=EvalRestriction(P6);assert Nullspace(AS)eq Nullspace(A6)and Rank(BasisMatrix(U)*AS)eq 6;assert Rank(BasisMatrix(U)*HorizontalJoin(AS,ColumnSubmatrix(AC,1,4)))eq 10;
print "C74_PARENT_COMPLETION_CUSP_AND_SEXTIC_GEOMETRY_PASS",18,Cputime(tt);keys:=[* *];
for i in[1..99]do row:=decs[i];P:=Prime(row[1]);assert Degree(P)eq 10;AA,rr:=EvalRestriction(P);assert Nullspace(AA)eq Nullspace(DM(row[2]))and Rank(BasisMatrix(U)*AA)eq 10 and row[3]eq 10;Append(~keys,Key2([Evaluate(xx,P),Evaluate(yy,P)]));if i mod 10 eq 0 then print "C74_PARENT_NATIVE_NC10_PROGRESS",i,Cputime(tt);end if;end for;
assert #Seqset([Sprint(v):v in keys])eq 99;Write(root cat"parent_native_point_keys.m",Sprint(keys):Overwrite:=true);Write(root cat"parent_native_summary_v2.m",Sprint([*99,18,10,38*]):Overwrite:=true);print "D31B_C74_PARENT_NATIVE_COMPLETE",99,18,Cputime(tt);quit;
