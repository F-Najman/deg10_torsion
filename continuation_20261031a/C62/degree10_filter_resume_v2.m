load "work/deg10/continuation_20261031a/C62/class_common_v2.m";
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
Write("work/deg10/continuation_20261031a/C62/degree10_resume_v2_census.m",Sprintf("N:=%o; PointKeys:=%o; EligibleTwoCounts:=%o; Representatives:=%o; Orbits:=%o;",N,keys,eligible,reps,orbits));
load "work/deg10/continuation_20261031a/C62/degree10_prefix_frozen.m";
function DE(e)return F![DecRat35(a):a in e];end function;
tests:=ClassTests; survivors:=[* *]; completed:={Integers()|};
for a in tests do
 i:=a[1];assert i notin completed;Include(~completed,i);Q:=qq[i];
 f1:=DE(a[4][1]);f2:=DE(a[4][2]);assert ideal<OF|f1,f2> eq Ideal(Q);
 extra,ffdata,ii,es:=Explode(a[5]);ff:=DE(ffdata);A:=(ExpectedGenus+10+extra)*Divisor(fb[ic0]);RR:=&+[Divisor(fb[ii[z]])*es[z]:z in [1..#ii]];
 assert Divisor(ff)+A-Divisor(Q) eq RR;
 cl:=[* *];for j in [1..#Mods] do Append(~cl,-&+[CG[j]|es[z]*CC[j][ii[z]]:z in [1..#ii]]);end for;
 assert [Eltseq(v):v in cl] eq a[2];yes:=&and[cl[j] in CH[j]:j in [1..#Mods]];assert yes eq a[3];if yes then Append(~survivors,a);end if;
end for;
Write("work/deg10/continuation_20261031a/C62/degree10_resume_v2_classes.m",Sprintf("ClassTests:=%o;",tests));
printf "D31A_C62_PREFIX_REPLAY_COMPLETE cases=%o time=%o\n",#tests,Cputime(tt);

for z in [1..#reps] do
 i:=reps[z];if eligible[i] eq 0 or i in completed then continue;end if;
 cl,wit:=SmoothClass(qq[i]);yes:=&and[cl[j] in CH[j]:j in [1..#Mods]];
 f1,f2:=TwoGenerators(qq[i]);row:=[*i,[Eltseq(c):c in cl],yes,[*FER(f1),FER(f2)*],wit*];Append(~tests,row);if yes then Append(~survivors,row);end if;
 Write("work/deg10/continuation_20261031a/C62/degree10_resume_v2_classes.m",Sprintf("// Case %o\nAppend(~ClassTests,%o);",i,row));
 printf "D31A_NEW_CLASS level=%o representative=%o orbit=%o/%o rational_image=%o time=%o\n",N,i,z,#reps,yes,Cputime(tt);
end for;
Write("work/deg10/continuation_20261031a/C62/degree10_resume_v2_complete.m",Sprintf("N:=%o; Tests:=%o; Survivors:=%o;",N,tests,[x[1]:x in survivors]));
printf "D31A_DEGREE10_CLASS_FILTER_COMPLETE level=%o points=%o tests=%o survivors=%o time=%o\n",N,#qq,#tests,#survivors,Cputime(tt);quit;
