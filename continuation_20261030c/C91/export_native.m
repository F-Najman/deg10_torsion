load "work/deg10/continuation_20261030c/C91/tate7_quotient_common.m";
function RatER(q)
 return [Eltseq(Numerator(q)),Eltseq(Denominator(q))];
end function;
function FER(q)
 return [RatER(FB!Coefficient(q,i)):i in [0..13]];
end function;
function DecodeRat(e)
 return Evaluate(PolynomialRing(GF(2))!e[1],t)/Evaluate(PolynomialRing(GF(2))!e[2],t);
end function;
function Decode(e)
 return &+[CY!DecodeRat(e[i])*s^(i-1):i in [1..14]];
end function;
qs:=[Q:Q in Zeros(CY!(t^2+t+1))|Degree(Q) eq 2];assert #qs eq 2;
allp:=cs cat qs; primeER:=[];
for i in [1..#allp] do
 P:=allp[i];ff,gg:=TwoGenerators(P);aa:=FER(ff);bb:=FER(gg);assert Decode(aa) eq ff and Decode(bb) eq gg;
 assert Valuation(ff,P) gt 0 and Valuation(gg,P) gt 0;
 Append(~primeER,[*Degree(P),aa,bb*]);
 if i gt #cs then
  u:=UniformizingElement(P); lam:=Evaluate(tau(u)/u,P);assert lam^3 eq 1 and lam ne 1;
  printf "C91_CM_TAME_DERIVATIVE CMindex=%o lambda_order3=true\n",i-#cs;
 end if;
end for;
Pfile:="work/deg10/continuation_20261030c/C91/native_packet.m";
Write(Pfile,Sprintf("DefPolyER:=%o; PrimeER:=%o;",[RatER(c):c in Coefficients(f)],primeER));
VO:=VectorSpace(GF(2),23); records:=[];
for i in [1..#cs] do for m in [1..10 div Degree(cs[i])] do
 V,vm:=DifferentialSpace(m*Divisor(cs[i]));W:=sub<VO|[VO!Eltseq(vm(v)@@om):v in Basis(V)]>;
 Append(~records,<1,i,m,0,0,Dimension(W),Eltseq(BasisMatrix(W))>);
end for;end for;
ratidx:=[i:i in [1..#cs]|Degree(cs[i]) eq 1];
for j in [1..#qs] do for a in [1..3] do for b in [a..3] do
 C:=4*Divisor(qs[j])+Divisor(cs[ratidx[a]])+Divisor(cs[ratidx[b]]);
 V,vm:=DifferentialSpace(C);W:=sub<VO|[VO!Eltseq(vm(v)@@om):v in Basis(V)]>;
 assert Dimension(U)-Dimension(U meet W) eq 8;
 Append(~records,<2,#cs+j,4,ratidx[a],ratidx[b],Dimension(W),Eltseq(BasisMatrix(W))>);
end for;end for;end for;
Write(Pfile,Sprintf("KernelRecords:=%o; CM:=Matrix(GF(2),23,23,%o); MT:=Matrix(GF(2),23,23,%o); UB:=Matrix(GF(2),%o,23,%o);",records,Eltseq(CM),Eltseq(MT),Dimension(U),Eltseq(BasisMatrix(U))));
printf "C91_NATIVE_PACKET_EXPORT_COMPLETE primes=%o kernels=%o time=%o\n",#allp,#records,Cputime(tt);quit;
