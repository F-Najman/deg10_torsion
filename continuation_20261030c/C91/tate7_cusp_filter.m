load "work/deg10/continuation_20261030c/C91/tate7_quotient_common.m";

assert Dimension(U) ge 10;
VO:=VectorSpace(F,23); localjets:=AssociativeArray(); kernels:=AssociativeArray();
outfile:="work/deg10/continuation_20261030c/C91/tate7_cusp_blocks.m";
Write(outfile,Sprintf("GenusY:=23; ds:=%o; UB:=Matrix(GF(2),%o,23,%o); localjets:=AssociativeArray(); kernelspaces:=AssociativeArray();",[Degree(P):P in cs],Dimension(U),Eltseq(BasisMatrix(U))));
for i in [1..#cs] do
 for m in [1..10 div Degree(cs[i])] do
  V,vm:=DifferentialSpace(m*Divisor(cs[i]));
  W:=sub<VO|[VO!Eltseq(vm(v)@@om):v in Basis(V)]>;
  A:=Transpose(BasisMatrix(Nullspace(Transpose(BasisMatrix(W)))));
  assert Ncols(A) eq m*Degree(cs[i]);
  localjets[<i,m>]:=BasisMatrix(U)*A;
  Write(outfile,Sprintf("kernelspaces[<%o,%o>]:=Matrix(GF(2),%o,23,%o); localjets[<%o,%o>]:=Matrix(GF(2),%o,%o,%o);",i,m,Dimension(W),Eltseq(BasisMatrix(W)),i,m,Nrows(localjets[<i,m>]),Ncols(localjets[<i,m>]),Eltseq(localjets[<i,m>])));
  printf "C91_TATE7_CUSP_BLOCK cusp=%o multiplicity=%o rank=%o time=%o\n",i,m,Rank(localjets[<i,m>]),Cputime(tt);
 end for;
end for;
function Weighted(ds,d)
 if #ds eq 0 then return d eq 0 select [[]] else []; end if;
 return &cat[[[k] cat v:v in Weighted(ds[2..#ds],d-k*ds[1])]:k in [0..d div ds[1]]];
end function;
vs:=Weighted([Degree(P):P in cs],10); assert #vs eq 120;
bad:=[];
for v in vs do
 A:=ZeroMatrix(F,Dimension(U),0);
 for i in [1..#cs] do if v[i] gt 0 then A:=HorizontalJoin(A,localjets[<i,v[i]>]); end if;end for;
 assert Ncols(A) eq 10;
 if Rank(A) ne 10 then Append(~bad,v);printf "C91_TATE7_CUSP_EXCEPTION vector=%o rank=%o\n",v,Rank(A);end if;
end for;
Write(outfile,Sprintf("tested:=%o; bad:=%o;",#vs,bad));
printf "C91_TATE7_CUSP_FILTER_COMPLETE patterns=%o exceptions=%o time=%o\n",#vs,#bad,Cputime(tt);
quit;
