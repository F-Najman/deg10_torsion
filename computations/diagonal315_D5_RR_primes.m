/* Necessary three-multiple diagonal test. ALL F4 effective degree4 divisors
   are included. AtF7 the entirely cuspidal branch has a global F-cusp lift
   and is excluded separately using odd-prime injectivity and gon_F>=6.
   Noncuspidal F7 patterns are single2/3/4 with every cusp complement,
   and2+2 including repetition. No global full-generation claim is used. */
SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(4*10^9);t0:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
Attach("work/deg10/shared/deps/class_rational_reduction.m");
Free:=FreeAbelianGroup(20);Xs:=[* *];Gs:=[* *];Hs:=[* *];hs:=[* *];ms:=[* *];
css:=[* *];basecs:=[* *];signatures:=[];degrees:=[];
function WeightedVectors(ds,d)
 if #ds eq 0 then return d eq 0 select [[]] else [];end if;
 result:=[];
 for k in [0..d div ds[1]] do
  for v in WeightedVectors(ds[2..#ds],d-k*ds[1]) do Append(~result,[k] cat v);end for;
 end for;
 return result;
end function;
lift:=eval Read("work/deg10/exclusions/315d5/unknown_cuspidal_lifts.txt");results:=[* *];
for p in [13,19,31,37] do
 X:=MDX11(3,15,GF(p));C:=Curve(X);FF:=FunctionField(C);
 E:=EllipticCurve([FF!f:f in X`_E]);P:=E![FF!f:f in X`_P];Q:=E![FF!f:f in X`_Q];
 ai:=aInvariants(E);assert ai[2] eq 0 and ai[4] eq 0 and ai[5] eq 0;
 assert P eq E![0,0];fs:=[ai[3],ai[1]^3-27*ai[3]];
 for k in [1..7] do
  for j in [0..2] do
   A:=k*Q+j*P;assert A[3] ne 0 and A[1] ne 0 and A[2] ne 0;
   Append(~fs,A[1]);Append(~fs,A[2]);
  end for;
  A:=k*Q;B:=A+P;assert A ne B and A ne -B;Append(~fs,A[1]-B[1]);
 end for;
 groups:=AssociativeArray();cs:=Cusps(X);
 for cusp in cs do
  sig:=[Valuation(f,cusp):f in fs];
  if IsDefined(groups,sig) then Append(~groups[sig],cusp);else groups[sig]:=[cusp];end if;
 end for;
 sigs:=Sort(Setseq(Keys(groups)));assert #sigs eq 20;
 blocks:=[&+[Divisor(c):c in groups[sig]]:sig in sigs];
 ds:=[Degree(o):o in blocks];assert {*d:d in ds*} eq {*1^^16,4^^4*};
 assert sigs eq lift[1] and ds eq lift[2];
 c0:=Support(blocks[Index(ds,1)])[1];
 print "DECIC315_DIAGONAL_CUSP_BLOCKS",p,#sigs,ds,Cputime(t0);

 dims:=[];patterns:=[* *];
 for k in [1..#lift[3]] do
  v:=lift[3][k];D:=5*Divisor(c0)+&+[v[i]*(blocks[i]-ds[i]*Divisor(c0)):i in [1..20]];assert Degree(D) eq 5;
  V,vm:=RiemannRochSpace(D);Append(~dims,Dimension(V));
 end for;
 Append(~results,[*p,dims*]);print "DECIC315_D5_CUSP_CLASS_EFFECTIVITY",p,[#{ k:k in [1..#dims]|dims[k] eq j}:j in [0..Maximum(dims)]],Cputime(t0);
end for;
Write("work/deg10/exclusions/315d5/cuspidal_class_RR_checks.txt",Sprint(results));
left:=[k:k in [1..#lift[3]]|&and[r[2][k] gt 0:r in results]];
print "DECIC315_D5_CUSP_CLASS_PRIME_SURVIVORS",left;print "DECIC315_D5_CUSP_CLASS_RR_COMPLETE";quit;
