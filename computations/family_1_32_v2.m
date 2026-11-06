SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(4*10^9);t0:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
X:=MDX1(32,Rationals());C:=Curve(X);F:=BaseRing(C);d:=10;
printf "DECIC_MODEL m=1 n=32 genus=%o base_degree=%o target_degree=%o time=%o\n",Genus(C),Degree(F),d,Cputime(t0);
if Genus(C) lt d then
 printf "DECIC_GENUS_FAMILY m=1 n=32 relative_degree=%o absolute_degree=10 rational_cusp=true genus=%o\n",d,Genus(C);
 print "DECIC_GENUS_FAMILY_COMPLETE";quit;
end if;
FF:=FunctionField(C);
try
 uu:=FF!(CoordinateRing(C).1);vv:=FF!(CoordinateRing(C).2);
 trial:=[uu,vv,uu+vv,uu-vv,uu*vv,uu/vv,uu^2/vv,vv^2/uu];
 for ff in trial do
  if Degree(PoleDivisor(ff)) eq d then
   printf "DECIC_PROJECTION_FAMILY m=1 n=32 relative_degree=%o absolute_degree=10 function=%o\n",d,ff;
   print "DECIC_PROJECTION_FAMILY_COMPLETE";quit;
  end if;
 end for;
catch err
 print "DECIC_PROJECTION_TRIAL_UNAVAILABLE";
end try;
os:=[Type(O) eq PlcCrvElt select Divisor(O) else O:O in CuspOrbitsQ(X)];
cs:=[O:O in os|Degree(O) eq 1];printf "DECIC_MODEL_CUSPS m=1 n=32 rational=%o\n",#cs;
found:=false;
for i in [1..Minimum(6,#cs)] do
 V,vm:=RiemannRochSpace(d*cs[i]);
 for v in Basis(V) do
  f:=vm(v);if IsZero(f) then continue;end if;
  if Degree(PoleDivisor(f)) eq d then
   printf "DECIC_FUNCTION_FAMILY m=1 n=32 relative_degree=%o absolute_degree=10 cusp_index=%o function=%o\n",d,i,f;found:=true;break;
  end if;
 end for;
 if found then break;end if;
end for;
printf "DECIC_MODEL_PROBE_COMPLETE m=1 n=32 found=%o time=%o\n",found,Cputime(t0);quit;
