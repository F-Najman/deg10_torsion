SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(4*10^9);tt:=Cputime();load "work/deg10/C236_resolution_20261030/cuspidal_boundary/lattice.m";
L:=LatticeWithBasis(Matrix(Integers(),CuspUnitBasis));L:=LLL(L);assert Rank(L) eq 35;print "C236_CUSP_LLL_READY",Minimum([Norm(v):v in Basis(L)]),Cputime(tt);
vs:=ShortVectors(L,104);out:=[* *];good:=[* *];
for r in vs do
 v:=Eltseq(r[1]);norm:=r[2];assert norm le 104;pos:=&+[x:x in v|x gt 0];neg:=-&+[x:x in v|x lt 0];assert pos eq neg;
 if pos eq 10 then
  counts:=[#[o:o in CuspOrbits|#o eq 1 and v[o[1]+1] gt 0],#[o:o in CuspOrbits|#o eq 1 and v[o[1]+1] lt 0]];Sort(~counts);n1:=counts[2];n2:=counts[1];hit:=n1 ge 4 and n1+5*n2 ge 24;Append(~out,[*v,norm,counts,hit*]);if hit then Append(~good,[*v,norm,counts*]);end if;
 end if;
end for;
Write("work/deg10/C236_resolution_20261030/cuspidal_boundary/short_vectors.m",Sprint([*#vs,out,good*]):Overwrite:=true);print "C236_CUSP_SHORT_COMPLETE",#vs,#out,#good,Cputime(tt);quit;
