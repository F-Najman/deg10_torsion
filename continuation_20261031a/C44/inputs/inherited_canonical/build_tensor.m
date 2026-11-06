SetColumns(0);SetQuitOnError(true);SetMemoryLimit(20*10^9);tt:=Cputime();
load "work/degree9/cyclic/continuation_20261026/C44_canonical_ring/canonical_data.m";
P:=ProjectiveSpace(GF(3),35);R:=CoordinateRing(P);mons:=[R.(p[1])*R.(p[2]):p in pairs];eqs:=[&+[v[i]*mons[i]:i in [1..666]]:v in Rows(QuadRelations)];IC:=ideal<R|eqs>;GB:=GroebnerBasis(IC);LM:=[LeadingMonomial(f):f in GB];
function StdMons(dd)
 ss:={R!1};for k in [1..dd] do cand:={m*R.i:m in ss,i in [1..36]};ss:={m:m in cand|not exists{l:l in LM|IsDivisibleBy(m,l)}};assert #ss eq (k eq 1 select 36 else (2*k-1)*35);end for;return Sort(Setseq(ss));
end function;
S2:=StdMons(2);S4:=StdMons(4);S6:=StdMons(6);assert #S2 eq 105 and #S4 eq 245 and #S6 eq 385;SS6:={@m:m in S6@};
out:=Open("work/degree9/cyclic/continuation_20261026/C44_canonical_ring/multiplication_data.m","w");fprintf out,"Std2Exponents:=%o;Std4Exponents:=%o;Std6Exponents:=%o;TensorBlocks:=[* *];\n",[Exponents(m):m in S2],[Exponents(m):m in S4],[Exponents(m):m in S6];Flush(out);
for i in [1..105] do rows:=[];for j in [1..245] do nf:=NormalForm(S2[i]*S4[j],IC);rr:=[GF(3)!0:k in [1..385]];mm:=Monomials(nf);cc:=Coefficients(nf);for k in [1..#mm] do ii:=Index(SS6,mm[k]);assert ii gt 0;rr[ii]:=cc[k];end for;Append(~rows,rr);end for;MM:=Matrix(GF(3),rows);fprintf out,"Append(~TensorBlocks,Matrix(GF(3),245,385,%o));\n",[Integers()!x:x in Eltseq(MM)];Flush(out);if i mod 10 eq 0 then printf "NONIC44_RING_TENSOR_BLOCK index=%o time=%o\n",i,Cputime(tt);end if;end for;
fprintf out,"assert #TensorBlocks eq 105;\n";delete out;printf "NONIC44_RING_MULTIPLICATION_CERTIFIED R2=105 R4=245 R6=385 products25725=true time=%o\n",Cputime(tt);quit;
