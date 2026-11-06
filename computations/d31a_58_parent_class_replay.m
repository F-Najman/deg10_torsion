SetColumns(0);SetNthreads(1);SetQuitOnError(true);Alarm(1400);tt:=Cputime();
load "work/deg10/continuation_20261031a/C58/class_common.m";
load "work/deg10/continuation_20261031a/C58/degree10_complete.m";
load "work/deg10/continuation_20261031a/C58/degree10_census.m";
assert #Tests eq 293;out:=[];surv:=[];
for row in Tests do
 idx,oldclasses,oldyes,ep,wit:=Explode(row);ff:=DecElt35(ep[1]);gg:=DecElt35(ep[2]);I:=ideal<OF|ff,gg>;assert IsPrime(I);Q:=Place(I);assert Degree(Q) eq 10 and Valuation(ff,Q) gt 0 and Valuation(gg,Q) gt 0;
 assert NewKey([Evaluate(v,Q):v in coords]) eq Sprint(PointKeys[idx]);
 extra,ef,ii,ee:=Explode(wit);a:=DecElt35(ef);assert a ne 0 and extra ge 0 and #ii eq #ee and &and[e gt 0:e in ee];assert &+[ee[j]*Degree(fb[ii[j]]):j in[1..#ii]] eq ExpectedGenus+extra;
 A:=(ExpectedGenus+10+extra)*Divisor(fb[ic0])-Divisor(Q);IA,JA:=Ideals(A);IR:=ideal<OF|1>;JR:=ideal<OI|1>;
 for j in[1..#ii] do assert ii[j] in[1..#fb];if IsFinite(fb[ii[j]])then IR*:=Ideal(fb[ii[j]])^ee[j];else JR*:=Ideal(fb[ii[j]])^ee[j];end if;end for;
 assert (a*OF)*IA eq IR and(a*OI)*JA eq JR;
 newclasses:=[-&+[CG[j]|ee[z]*CC[j][ii[z]]:z in[1..#ii]]:j in[1..#Mods]];
 assert &and[newclasses[j] eq CG[j]!oldclasses[j]:j in[1..#Mods]];yes:=&and[newclasses[j]in CH[j]:j in[1..#Mods]];assert yes eq oldyes;
 if yes then Append(~surv,idx);end if;Append(~out,[idx,extra,#ii,yes select 1 else 0]);
 if #out mod 25 eq 0 then print "C58_PARENT_CLASS_PROGRESS",#out,Cputime(tt);end if;
end for;
assert #surv eq 2 and Seqset(surv) eq {462,1748};Write("work/deg10/continuation_20261031a/C58/parent_class_replay.m",Sprint([*out,surv*]):Overwrite:=true);print "C58_PARENT_CLASS_REPLAY_COMPLETE",#out,surv,Cputime(tt);quit;
