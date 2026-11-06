load "work/deg10/continuation_20261031a/C62/class_common.m";
Write("work/deg10/continuation_20261031a/C62/degree10_classes.m","ClassTests:=[* *];");

AttachSpec("work/degree9/shared/mdmagma/v2/mdmagma.spec");
qq:=Places(F,10); expected:=5685;assert #qq eq expected;
keys:=[];PointIndex:=AssociativeArray();allXY:=[* *];eligible:=[];
for i in [1..#qq] do
 Q:=qq[i];assert Degree(Q) eq 10 and Q notin cs;
 xy:=[Evaluate(a,Q):a in coords];ky:=NewKey(xy);assert not IsDefined(PointIndex,ky);PointIndex[ky]:=i;Append(~keys,ky);Append(~allXY,xy);
 bv:=Evaluate(bb,Q);cv:=Evaluate(cc,Q);L:=Parent(bv);EL:=EllipticCurve([L|1-cv,-bv,-bv,0,0]);PP:=EL![0,0];assert N*PP eq EL!0 and PP ne EL!0;
 R<Z>:=PolynomialRing(L);ai:=aInvariants(EL);two:=4*Z^3+(ai[1]^2+4*ai[2])*Z^2+2*(ai[1]*ai[3]+2*ai[4])*Z+ai[3]^2+4*ai[5];nr:=#Roots(two);assert nr in [0,1,3];Append(~eligible,nr);
 if i mod 1000 eq 0 then printf "D31A_NEW_DEG10_POINTS_READY level=%o checked=%o/%o time=%o\n",N,i,#qq,Cputime(tt);end if;
end for;
printf "D31A_NEW_DEG10_POINTS_COMPLETE level=%o points=%o rational2count=%o time=%o\n",N,#qq,#[a:a in eligible|a gt 0],Cputime(tt);
visited:={};reps:=[];orbits:=[* *];
for i in [1..#qq] do
 if i in visited then continue;end if;
 Q:=qq[i];bv:=Evaluate(bb,Q);cv:=Evaluate(cc,Q);L:=Parent(bv);EL:=EllipticCurve([L|1-cv,-bv,-bv,0,0]);PP:=EL![0,0];orb:={};
 for a in [1..(N-1) div 2] do
  ky:=NewKey(TateNormalForm_xy(EL,a*PP));assert IsDefined(PointIndex,ky);Include(~orb,PointIndex[ky]);
 end for;
 assert #orb gt 0 and (N-1) mod (2*#orb) eq 0 and &and[j notin visited and eligible[j] eq eligible[i]:j in orb];
 visited join:=orb;Append(~reps,i);Append(~orbits,Sort(Setseq(orb)));
end for;assert #visited eq #qq;
printf "D31A_NEW_DIAMOND_ORBITS level=%o orbits=%o eligible_orbits=%o time=%o\n",N,#reps,#[i:i in reps|eligible[i] gt 0],Cputime(tt);
Write("work/deg10/continuation_20261031a/C62/degree10_census.m",Sprintf("N:=%o; PointKeys:=%o; EligibleTwoCounts:=%o; Representatives:=%o; Orbits:=%o;",N,keys,eligible,reps,orbits));
tests:=[* *];survivors:=[* *];
for z in [1..#reps] do
 i:=reps[z];if eligible[i] eq 0 then continue;end if;
 cl,wit:=SmoothClass(qq[i]);yes:=&and[cl[j] in CH[j]:j in [1..#Mods]];
 f1,f2:=TwoGenerators(qq[i]);row:=[*i,[Eltseq(c):c in cl],yes,[*FER(f1),FER(f2)*],wit*];Append(~tests,row);if yes then Append(~survivors,row);end if;
 Write("work/deg10/continuation_20261031a/C62/degree10_classes.m",Sprintf("// Case %o\nAppend(~ClassTests,%o);",i,row));
 printf "D31A_NEW_CLASS level=%o representative=%o orbit=%o/%o rational_image=%o time=%o\n",N,i,z,#reps,yes,Cputime(tt);
end for;
Write("work/deg10/continuation_20261031a/C62/degree10_complete.m",Sprintf("N:=%o; Tests:=%o; Survivors:=%o;",N,tests,[x[1]:x in survivors]));
printf "D31A_DEGREE10_CLASS_FILTER_COMPLETE level=%o points=%o tests=%o survivors=%o time=%o\n",N,#qq,#tests,#survivors,Cputime(tt);quit;
