SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(2*10^9);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");Attach("work/deg10/continuation_20261029b/deps/certified_kernel_isogenies.m");
data:=eval Read("work/deg10/continuation_20261029b/C250/filter25.m");out:=[* *];
function Code(v) return [Integers()!c:c in Eltseq(v,GF(3))];end function;
for desc in data[7] do
 d:=desc[2];row:=[r:r in data[2]|r[1] eq d][desc[3]];ky:=row[3];a:=[r:r in data[4]|r[1] eq d][1];R:=PolynomialRing(GF(3));K<z>:=ext<GF(3)|R!a[2]>;b:=Evaluate(R!ky[1..d],z);c:=Evaluate(R!ky[d+1..2*d],z);E:=EllipticCurve([K|0,c,0,b*(1-b-c),0]);Q:=E![b,b];assert Order(Q) eq 50;
 isogs:=MDCertifiedIsogenies(E,Q,7);rows:=[* *];
 for rec in isogs do
  EG,PG,e,g,U,kg,Ng:=Explode(rec);L:=BaseRing(EG);Append(~rows,[*e,[Integers()!c:c in Coefficients(DefiningPolynomial(L,GF(3)))],Code(kg),[Code(v):v in Coefficients(g)],[Code(v):v in aInvariants(EG)],[Code(PG[1]),Code(PG[2])],[Code(v):v in Coefficients(Ng)]*]);
 end for;
 assert &+[r[1]:r in rows] eq 8;Append(~out,[*desc,d,a[2],ky,rows*]);print "C250_LOCAL_KERNEL_ORBITS_COMPLETE",desc,[r[1]:r in rows],Cputime(tt);
end for;
Write("work/deg10/continuation_20261029b/C250/kernel_certificates.m",Sprint(out));print "ALL_C250_LOCAL_KERNEL_CERTIFICATES_COMPLETE",Cputime(tt);quit;
