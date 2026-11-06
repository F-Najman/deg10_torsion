SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(8*10^9);Alarm(480);tt:=Cputime();AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");old:=eval Read("work/degree8/noncyclic/N28_D8_quotient_cusp_RR_arithmetic.txt");sigs,ds,base,lifts:=Explode(old[1..4]);out:=[* *];
for p in[5,11,13] do
 Y:=MDX1(28,GF(p));raw:=CuspOrbitsQ(Y);ss:=[CuspSignature(Y,Support(o)[1]):o in raw];ord:=Sort([<ss[i],i>:i in[1..#raw]]);os:=[raw[z[2]]:z in ord];assert [ss[z[2]]:z in ord] eq sigs and [Degree(o):o in os] eq ds;c0:=Support(os[base])[1];L:=8*Divisor(c0)+&+[lifts[1][i]*(os[i]-ds[i]*Divisor(c0)):i in[1..15]];assert Degree(L) eq 8;rats:=[o:o in os|Degree(o) eq 1];quads:=[o:o in os|Degree(o) eq 2];cs:=[rats[i]+rats[j]:i,j in[1..#rats]|i le j] cat quads;assert #cs eq 48;hs:=[];
 for i in[1..#cs] do D:=FunctionFieldDivisor(L+cs[i]);Append(~hs,#Basis(D));end for;
 Append(~out,[*p,hs*]);print "Y228_OCTIC8_L_PLUS_C2_RR",p,hs,Cputime(tt);Write("work/deg10/continuation_20261031b/Y228_shared/octic8/L8_plus_cusp2_RR.m",Sprint(out):Overwrite:=true);
 if &and[h eq 1:h in hs] then break;end if;
end for;print "D31B_Y228_OCTIC8_SMALL_RR_COMPLETE",Cputime(tt);quit;
