SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(5*10^9);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
a:=eval Read("work/deg10/continuation_20261029b/C226/small_cusp_divisors_v2.m");
X:=MDX11(2,26,Rationals());C:=Curve(X);assert Genus(C) eq 25;print "C226_Q_DIRECT_MODEL_READY",Cputime(tt);b:=X`_coordinates[1];pp,ee:=Support(Divisor(b));cs:=[c:c in pp|Degree(c) eq 1];labels:=[Valuation(b,c):c in cs];assert #cs eq 18 and Seqset(labels) eq {v:v in a[1]};out:=[* *];
for row in a[2] do
 DD:=&+[row[3][i]*Divisor(cs[Index(labels,a[1][i])]):i in [1..18]];assert Degree(DD) eq 10;print "C226_Q_DIRECT_RR_START",row[1],Cputime(tt);
 V,mu:=RiemannRochSpace(DD);dim:=Dimension(V);print "C226_Q_DIRECT_RR",row[1],dim,Cputime(tt);Append(~out,[*row[1],row[2],dim*]);
 if dim eq 1 then hh:=mu(V.1);EE:=DD+Divisor(hh);assert IsEffective(EE) and Degree(EE) eq 10;ps,ms:=Support(EE);print "C226_Q_DIRECT_EFFECTIVE_SUPPORT",row[1],[Degree(p):p in ps],ms,[p in cs:p in ps];Write(Sprintf("work/deg10/continuation_20261029b/C226/Q_direct_effective_%o.txt",row[1]),Sprint(EE) cat "\n" cat Sprint(hh));end if;
 Write("work/deg10/continuation_20261029b/C226/Q_direct_RR.m",Sprint(out):Overwrite:=true);
end for;
print "C226_Q_DIRECT_RR_COMPLETE",Cputime(tt);quit;
