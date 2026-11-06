SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(5*10^9);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
Attach("work/deg10/continuation_20261029/deps/fast_hecke.m");
Attach("work/deg10/continuation_20261029/deps/fast_hecke_11.m");
Attach("work/deg10/continuation_20261029/deps/fast_diamond_hash_m2.m");
Attach("work/deg10/continuation_20261029/deps/fast_full_aut_hash_m4_16.m");
Attach("work/deg10/continuation_20261029/deps/fast_full_aut_hash_m4_16_relative.m");
Attach("work/deg10/continuation_20261029/deps/full_aut_4_20.m");
X:=MDX11(4,20,GF(9));xs:=NoncuspidalPlaces(X,2);ys:=NoncuspidalPlaces(X,3);rr,labels,keys:=Decic420Orbits(X,xs);todo:=[];
for sh in [0..7] do
 old:=eval Read(Sprintf("work/deg10/continuation_20261029/cont_sieve_4_20_mixed_shard%o_arithmetic.txt",sh));assert keys eq old[1] and labels eq old[2] and [Index(xs,x):x in rr] eq old[3];
 for r in old[6] do assert r[1] eq 22;Append(~todo,[r[2],r[3]]);end for;
end for;assert #todo eq 8;tests:=[];left:=[];cache:=AssociativeArray();
for pair in todo do
 i,j:=Explode(pair);a:=Index(xs,rr[i]);killed:=false;
 for q in [13,17] do
  for k in Seqset([a,j]) do
   ky:=[q,k];if not IsDefined(cache,ky) then cache[ky]:=HeckeOperatorFast(X,q,xs[k])-Divisor(xs[k])-q*Divisor(DiamondOperator(X,q,xs[k]));end if;
  end for;
  pr:=IsPrincipal(cache[[q,a]]+cache[[q,j]]);Append(~tests,[i,j,q,pr select 1 else 0]);print "DECIC420_PAIR_FOLLOWUP",i,j,q,pr,Cputime(tt);
  if not pr then killed:=true;break;end if;
 end for;if not killed then Append(~left,pair);end if;
end for;
Write("work/deg10/continuation_20261029/M4N20_pair_followup_arithmetic.txt",Sprint([*keys,labels,[Index(xs,x):x in rr],tests,left*]));
function Images(X,x)
 F:=BaseRing(X);dd:=Degree(x);K:=ext<F|dd>;q,t:=Explode(M4SixteenRelativeToCommon([Evaluate(f,x):f in X`_coordinates],K,F));
 B:=(q^2-1)*(t^2-1)/16;E:=EllipticCurve([K!1,B,B,K!0,K!0]);P:=E![0,0];Q:=E![(q+1)*(t^2-1)/8,(q+1)^2*(t-1)^2*(t+1)/32];
 out:=[];
 for par in Decic420Parameters() do a,b,c,d:=Explode(par);Append(~out,M4SixteenRelativeKey(Decic420Normalize(E,a*P+5*b*Q,c*P+d*Q),K,F,dd));end for;return out;
end function;
pars:=Decic420Parameters();stabs:=[* *];
procedure AddStabilizer(typ,descr,ks,ims,~stabs)
 stabilizer:=[g:g in [1..96]|ims[g] eq ks];free:=[];
 for g in stabilizer do
  a,b,c,d:=Explode(pars[g]);u:=d mod 5;A:=Matrix(Integers(4),2,2,[a,b,c,d]);I:=IdentityMatrix(Integers(4),2);
  if u in [1,4] and A^2 eq I and not ((u eq 1 and A eq I) or (u eq 4 and A eq -I)) then Append(~free,g);end if;
 end for;assert #free gt 0;
 Append(~stabs,[*typ,descr,ks,ims,stabilizer,free*]);print "DECIC420_CUSPIDAL_CORRECTION_STABILIZER",typ,descr,#stabilizer,free,Cputime(tt);
end procedure;
for pair in left do
 i,j:=Explode(pair);aa:=Images(X,rr[i]);bb:=Images(X,xs[j]);ks:=Sort([keys[Index(xs,rr[i])],keys[j]]);ims:=[Sort([aa[g],bb[g]]):g in [1..96]];AddStabilizer(22,pair,ks,ims,~stabs);
end for;
xs4:=NoncuspidalPlaces(X,4);rs4,ll4,kk4:=Decic420Orbits(X,xs4);old:=eval Read("work/deg10/continuation_20261029/M4N20_quartic_followup_arithmetic.txt");assert kk4 eq old[2] and ll4 eq old[3];
for i in old[6] do
 im:=Images(X,rs4[i]);ks:=[kk4[Index(xs4,rs4[i])]];ims:=[[key]:key in im];AddStabilizer(4,[i],ks,ims,~stabs);
end for;
Write("work/deg10/continuation_20261029/M4N20_D5_stabilizers_arithmetic.txt",Sprint([*pars,stabs*]));print "DECIC420_PAIR_FOLLOWUP_AND_STABILIZERS_COMPLETE",left,#stabs,Cputime(tt);quit;
