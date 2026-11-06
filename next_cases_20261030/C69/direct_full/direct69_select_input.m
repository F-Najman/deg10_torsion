SetColumns(0);SetNthreads(1);SetQuitOnError(true);SetMemoryLimit(4*10^9);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
data:=eval Read("work/deg10/next_cases_20261030/C69/degree10_base.m");assert data[1] eq 20 and data[2] eq 2;
jcode:=eval Read("work/deg10/next_cases_20261030/C69/base_j_functions.m");assert #jcode eq 2;
F<a>:=GF(2^10);
function FieldValue(bits)
 return &+[F!bits[i]*a^(i-1):i in [1..#bits]];
end function;
function EvalRF(code,x)
 R<T>:=PolynomialRing(F);return Evaluate(R!code[1],x)/Evaluate(R!code[2],x);
end function;
CQ:=SmallModularCurve(23);jr:=jInvariant(CQ,23);CC:=ChangeRing(CQ,GF(2));RC:=CoordinateRing(Ambient(CC));jf:=FunctionField(CC)!((RC!Numerator(jr))/(RC!Denominator(jr)));
K,km:=AlgorithmicFunctionField(FunctionField(CC));assert Genus(K) eq 2;
numY,denY:=IntegralSplit(K.1,MaximalOrderFinite(K));assert numY eq K.1 and denY eq 1;
J:=km(jf);assert #Eltseq(J) eq 2;
S<T>:=PolynomialRing(GF(2));BX:=BaseField(K);
for i in [1,2] do assert BX!(S!jcode[i][1])/(BX!(S!jcode[i][2])) eq Eltseq(J)[i];end for;
jvalues:=[];jpolys:=[];
for idx in [1,2] do
 key:=data[3][idx][1];assert #key eq 30;
 xx:=FieldValue(key[1..10]);yy:=FieldValue(key[11..20]);
 kp:=DefiningPolynomial(K);assert &+[EvalRF([[Integers()!c:c in Coefficients(Numerator(v))],[Integers()!c:c in Coefficients(Denominator(v))]],xx)*yy^(i-1):i->v in Coefficients(kp)] eq 0;
 jj:=EvalRF(jcode[1],xx)+EvalRF(jcode[2],xx)*yy;Append(~jvalues,jj);Append(~jpolys,MinimalPolynomial(jj,GF(2)));assert Degree(jpolys[idx]) eq 10;
end for;
assert jpolys[1] ne jpolys[2];target:=jvalues[2];key:=data[3][2][1];choices:=[];
for E in Twists(EllipticCurveFromjInvariant(target)) do
 if #E mod 69 ne 0 then continue;end if;
 G,gm:=AbelianGroup(E);assert Invariants(G) eq [966];PN:=14*gm(G.1);assert Order(PN) eq 69;Append(~choices,<E,PN>);
end for;
assert #choices eq 1;E:=choices[1][1];PN:=choices[1][2];
fp:=[Integers()!c:c in Coefficients(DefiningPolynomial(F,GF(2)))];ecs:=[[Integers()!c:c in Eltseq(v,GF(2))]:v in aInvariants(E)];pcs:=[[Integers()!c:c in Eltseq(PN[k],GF(2))]:k in [1..2]];
record:=[*69,10,2,key,fp,ecs,pcs,[Integers()!c:c in Coefficients(jpolys[2])]*];
Write("work/deg10/next_cases_20261030/C69/direct_full/selected_input.m",Sprint(record):Overwrite:=true);
Write("work/deg10/next_cases_20261030/C69/direct_full/selected_hecke_input.m",Sprintf("N:=69;Ds:=[10];Qs:=[5];a:=1;Inputs:=[*<10,2,%o,%o,%o>*];",fp,ecs,pcs):Overwrite:=true);
printf "DIRECT69_SELECTED_INPUT_COMPLETE base_case=2 j_minpoly=%o cyclic_group=966 point_order=69 field_degree=10 exact_base_key=%o time=%o\n",jpolys[2],key,Cputime(tt);quit;
