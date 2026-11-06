SetColumns(0);SetQuitOnError(true);SetMemoryLimit(12*10^9);
AttachSpec("work/degree9/shared/mdmagma/v2/mdmagma.spec");
Attach("work/degree9/shared/fast_hecke.m");
N:=51;Ds:=[9,10];Qs:=[5];a:=1;RecordFile:="work/deg10/next_cases_20261030/cyclic_further/C51/finite_inputs.m";
t0:=Cputime();X:=MDX1(N,GF(2));C:=Curve(X);g:=Genus(C);assert g eq 65;
printf "DECIC51_MODEL N=%o genus=%o alpha_diamond=%o time=%o\n",N,g,a,Cputime(t0);
PrintFile(RecordFile,Sprintf("N:=%o;Ds:=%o;Qs:=%o;a:=%o;Inputs:=[**];",N,Ds,Qs,a):Overwrite:=true);
summary:=[];
for d in Ds do
 pl:=NoncuspidalPlaces(X,d);reps:=PlacesUpToDiamond(X,pl);
 printf "DECIC51_ENUM N=%o degree=%o all_places=%o diamond_reps=%o time=%o\n",N,d,#pl,#reps,Cputime(t0);
 survivors:=[];
 for i in [1..#reps] do
  P:=reps[i];D:=Divisor(P);E:=EllipticCurve(X,P);PN:=LevelStructure(X,P)`P;
  assert Order(PN) eq N;
  K:=BaseRing(E);assert Degree(K) eq d;
  PrintFile(RecordFile,Sprintf("Append(~Inputs,<%o,%o,%o,%o,%o>);",d,i,[Integers()!c:c in Eltseq(DefiningPolynomial(K,GF(2)))],[[Integers()!c:c in Eltseq(v,GF(2))]:v in aInvariants(E)],[[Integers()!c:c in Eltseq(PN[k],GF(2))]:k in [1..2]]));
  eliminated:=false;
  for q in Qs do
   assert IsPrime(q) and IsOdd(q) and N mod q ne 0;
   T:=HeckeOperatorFast(X,q,P);assert Degree(T) eq (q+1)*d;
   B:=T-D-q*Divisor(DiamondOperator(X,q,P));
   if a ne 1 then B:=DiamondOperator(X,a,B)-B;end if;
   assert Degree(B) eq 0;
   printf "DECIC51_TEST_START N=%o degree=%o index=%o q=%o time=%o\n",N,d,i,q,Cputime(t0);
   flag:=IsPrincipal(B);
   printf "DECIC51_TEST_RESULT N=%o degree=%o index=%o q=%o principal=%o time=%o\n",N,d,i,q,flag,Cputime(t0);
   if not flag then eliminated:=true;break;end if;
  end for;
  if not eliminated then Append(~survivors,i);end if;
 end for;
 printf "DECIC51_DEGREE_SUMMARY N=%o degree=%o all_places=%o diamond_reps=%o surviving_indices=%o time=%o\n",N,d,#pl,#reps,survivors,Cputime(t0);
 Append(~summary,<d,#pl,#reps,#survivors>);
end for;
printf "DECIC51_GOOD_BRANCHES_COMPLETED N=%o summary=%o time=%o\n",N,summary,Cputime(t0);quit;
