SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(8*10^9);Alarm(300);tt:=Cputime();root:="work/deg10/continuation_20261031b/C66/";AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
X:=MDX1(33,GF(5));FF:=FunctionField(Curve(X));AF,am:=AlgorithmicFunctionField(FF);E:=EllipticCurve([am(z):z in X`_E]);P:=E![AF|0,0];units:=[a:a in[1..16]|GCD(a,33) eq 1];os:=CuspOrbitsQ(X);blocks:=[FunctionFieldDivisor(o):o in os];sigs:=[CuspSignature(X,Support(o)[1]):o in os];assert #os eq 17;out:=[];
for a in units do
 g:=TateNormalForm_xy(E,a*P)[2];D:=PoleDivisor(g);assert Degree(D) eq 10;vs:=[];
 for B in blocks do ps,es:=Support(B);ww:=[Valuation(D,ps[i]) div es[i]:i in[1..#ps]];assert #Seqset(ww) eq 1;Append(~vs,ww[1]);end for;
 assert D eq &+[vs[i]*blocks[i]:i in[1..17]];Append(~out,vs);
end for;
Write(root cat "pencil_cusp_vectors_F5.m",Sprint([*sigs,units,out*]):Overwrite:=true);printf "D31B_C66_PENCIL_CUSP_CLASSES_COMPLETE all10degree10=true all_global_cusp_blocks=true time=%o\n",Cputime(tt);quit;
