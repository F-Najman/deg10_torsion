SetColumns(0);SetNthreads(1);SetQuitOnError(true);SetMemoryLimit(9*10^9);Alarm(650);tt:=Cputime();AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
root:="work/deg10/continuation_20261031a/C46/boundary/";X:=MDX1(46,GF(3));C:=Curve(X);assert Genus(C) eq 45;F,mp:=AlgorithmicFunctionField(FunctionField(C));assert Genus(F) eq 45;
cs:=Cusps(X);assert {*Degree(P):P in cs*} eq {*1^^22,11^^2*};rat:=[P:P in cs|Degree(P) eq 1];j:=mp(jInvariantMap(X));
ps:=[Support(DivisorGroup(F)!FunctionFieldDivisor(Divisor(P)))[1]:P in rat];widths:=[-Valuation(j,P):P in ps];assert {*w:w in widths*} eq {*23^^11,46^^11*};bases:=[Index(widths,23),Index(widths,46)];
V,om:=SpaceOfHolomorphicDifferentials(F);assert Dimension(V) eq 45;VO:=VectorSpace(GF(3),45);flags:=[];kernels:=[];
print "C46_BOUNDARY_MODEL_READY",Degree(F),Genus(F),bases,Cputime(tt);
for i in[1..#ps] do
 lim:=i in bases select 10 else 5;u:=UniformizingElement(ps[i]);du:=Differential(u);L,ev:=Completion(F,ps[i]:Precision:=160);series:=[ev(w/du):w in[om(v):v in Basis(V)]];
 assert &and[Valuation(w) ge 0 and AbsolutePrecision(w) ge lim:w in series];vv:=[[Integers()!Coefficient(w,n):w in series]:n in[0..lim-1]];
 for m in[1..lim] do A:=Matrix(GF(3),vv[1..m]);assert Rank(A) eq m;W:=Nullspace(Transpose(A));Append(~kernels,<i,m,Dimension(W),Eltseq(BasisMatrix(W))>);end for;
 Append(~flags,vv);print "C46_BOUNDARY_COMPLETION_FLAG",i,widths[i],lim,Cputime(tt);
 Write(root cat "completion_flags_progress.m",Sprint([*widths,bases,flags,kernels*]):Overwrite:=true);
end for;
Write(root cat "completion_restriction_flags.m",Sprint([*widths,bases,flags,kernels*]):Overwrite:=true);print "C46_BOUNDARY_COMPLETION_FLAGS_COMPLETE",#ps,#kernels,Cputime(tt);quit;
