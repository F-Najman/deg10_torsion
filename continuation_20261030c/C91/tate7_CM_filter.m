load "work/deg10/continuation_20261030c/C91/tate7_quotient_common.m";
qs:=[Q:Q in Zeros(CY!(t^2+t+1))|Degree(Q) eq 2]; assert #qs eq 2;
rat:=[P:P in cs|Degree(P) eq 1]; assert #rat eq 3;
VO:=VectorSpace(GF(2),23); tested:=0; bad:=[];
Write("work/deg10/continuation_20261030c/C91/tate7_CM_blocks.m", "CMRecords:=[];");
for j in [1..#qs] do
 f,g:=TwoGenerators(qs[j]);
 assert Valuation(tau(f),qs[j]) gt 0 and Valuation(tau(g),qs[j]) gt 0;
 for a in [1..3] do for b in [a..3] do
  Cref:=4*Divisor(qs[j])+Divisor(rat[a])+Divisor(rat[b]); assert Degree(Cref) eq 10;
  V,vm:=DifferentialSpace(Cref); W:=sub<VO|[VO!Eltseq(vm(v)@@om):v in Basis(V)]>;
  rr:=Dimension(U)-Dimension(U meet W); tested+:=1;
  Write("work/deg10/continuation_20261030c/C91/tate7_CM_blocks.m", Sprintf("Append(~CMRecords,<%o,%o,%o,%o,Matrix(GF(2),%o,23,%o)>);",j,a,b,rr,Dimension(W),Eltseq(BasisMatrix(W))));
  printf "C91_TATE7_CM_BLOCK CM=%o ratA=%o ratB=%o rank=%o expected8=true time=%o\n",j,a,b,rr,Cputime(tt);
  if rr ne 8 then Append(~bad,<j,a,b,rr>);end if;
 end for;end for;
end for;
assert tested eq 12;
printf "C91_TATE7_CM_FILTER_COMPLETE tested=%o exceptions=%o time=%o\n",tested,bad,Cputime(tt);quit;
