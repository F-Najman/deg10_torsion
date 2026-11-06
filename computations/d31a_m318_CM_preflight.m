SetColumns(0);SetNthreads(1);SetQuitOnError(true);SetMemoryLimit(8*10^9);Alarm(500);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");K<z>:=CyclotomicField(3);X:=MDX11(3,18,K);C:=Curve(X);assert Genus(C) eq 10;
F,mp:=AlgorithmicFunctionField(FunctionField(C));F:=RationalExtensionRepresentation(F);j:=F!mp(jInvariant(EllipticCurve(X`_E)));
print "M318_CM_NUMBER_FIELD_MODEL_READY",Genus(F),Degree(F),Cputime(tt);
ps:=Zeros(j);dd:=[Degree(P):P in ps];print "M318_CM_J0_DEGREES",dd,Cputime(tt);
Write("work/deg10/continuation_20261031a/M3N18/CM_preflight.m",Sprint([*dd,[Valuation(j,P):P in ps]*]):Overwrite:=true);print "M318_CM_NUMBER_FIELD_PREFLIGHT_COMPLETE",Cputime(tt);quit;
