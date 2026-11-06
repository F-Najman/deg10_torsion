SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(3*10^9);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
X:=MDX1(36,Rationals());C:=Curve(X);ff:=FunctionField(C);d:=10;u,v:=Explode(X`_coordinates_xy);
printf "DECIC_FUNCTION_START m=1 n=36 genus=%o target=10\n",Genus(C);
found:=false;
for i in [-4..4] do
 for j in [-4..4] do
  if i eq 0 and j eq 0 then continue;end if;
  f:=u^i*v^j;
  if Degree(f) eq d then
   printf "DECIC_FUNCTION_SUCCESS m=1 n=36 relative_degree=10 base_degree=%o exponents=%o function=%o\n",Degree(BaseRing(C)),[i,j],f;found:=true;break;
  end if;
 end for;
 if found then break;end if;
end for;
if not found then
 for f in [u+v,u-v,u+1/v,u-1/v,u^2+v,v^2+u,(u+1)/v,(v+1)/u] do
  if Degree(f) eq d then printf "DECIC_FUNCTION_SUCCESS m=1 n=36 relative_degree=10 base_degree=%o function=%o\n",Degree(BaseRing(C)),f;found:=true;break;end if;
 end for;
end if;
printf "DECIC_FUNCTION_SEARCH_COMPLETE m=1 n=36 found=%o time=%o\n",found,Cputime(tt);quit;
