SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(8*10^9);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");Attach("work/deg10/continuation_20261029/deps/model_6_12.m");
X:=MDX612(GF(13));C:=Curve(X);AF,af:=AlgorithmicFunctionField(FunctionField(C));print "DECIC612_AF_ENUM_BEGIN",Degree(AF),Genus(AF),Cputime(tt);
xs:=Places(AF,5);print "DECIC612_AF_ENUM_COUNT",#xs,Cputime(tt);assert #xs eq 74688;
Write("work/deg10/continuation_20261029/M6N12_AF13_enum_count.txt",Sprint([*Degree(AF),Genus(AF),#xs*]));print "DECIC612_AF_ENUM_COMPLETE",Cputime(tt);quit;
