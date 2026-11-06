SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(4*10^9);t0:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
X:=MDX1(2,12);C:=Curve(X);F:=BaseRing(C);d:=10;
printf "DECIC_MODEL m=2 n=24 genus=%o base_degree=%o target_degree=%o time=%o\n",Genus(C),Degree(F),d,Cputime(t0);
if Genus(C) lt d then
 printf "DECIC_GENUS_FAMILY m=2 n=24 relative_degree=%o absolute_degree=10 rational_cusp=true genus=%o\n",d,Genus(C);
end if;
os:=[Type(O) eq PlcCrvElt select Divisor(O) else O:O in CuspOrbitsQ(X)];
cs:=[O:O in os|Degree(O) eq 1];printf "DECIC_MODEL_CUSPS m=2 n=24 rational=%o\n",#cs;
found:=false;
for i in [1..Minimum(6,#cs)] do
 V,vm:=RiemannRochSpace(d*cs[i]);
 for v in Basis(V) do
  f:=vm(v);if IsZero(f) then continue;end if;
  if Degree(PoleDivisor(f)) eq d then
   printf "DECIC_FUNCTION_FAMILY m=2 n=24 relative_degree=%o absolute_degree=10 cusp_index=%o function=%o\n",d,i,f;found:=true;break;
  end if;
 end for;
 if found then break;end if;
end for;
printf "DECIC_MODEL_PROBE_COMPLETE m=2 n=24 found=%o time=%o\n",found,Cputime(t0);quit;
