SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(4*10^9);tt:=Cputime();root:="work/deg10/next_cases_20261030/C242/";
a:=eval Read(root cat "target14_proven_picard.m");Z:=FreeAbelianGroup(#a[6][1]);G,q:=quo<Z|[Z!r:r in a[6]]>;H:=sub<G|[q(Z!v):v in a[9]]>;assert Invariants(H) eq [2,2,6,18];Q,qh:=quo<G|H>;assert #TorsionSubgroup(Q) eq 4;
own:=AssociativeArray();small:=AssociativeArray();fields:=[* *];
for d in [4..9] do tag:=d le 6 select "small" else Sprintf("D%o",d);a:=eval Read(root cat Sprintf("filter14_D%o_%o.m",d,tag));Append(~fields,[*d,a[1]*]);own[d]:=a[2];if d le 6 then small[d]:=a[3];end if;end for;
dict:=AssociativeArray();field10:=[];
for s in [0..3] do a:=eval Read(root cat Sprintf("filter14_D10_D10s%o.m",s));if s eq 0 then field10:=a[1];else assert a[1] eq field10;end if;for r in a[2] do assert not IsDefined(dict,r[2]);dict[r[2]]:=r;end for;end for;
src10:=eval Read(root cat "source10.m");n10:=#[r:r in src10[9]|r[1] eq 10];assert Sort(Setseq(Keys(dict))) eq [1..n10];own[10]:=[*dict[i]:i in [1..n10]*];Append(~fields,[*10,field10*]);
assert [#own[d]:d in [4..9]] eq [5,18,64,316,1372,5888];assert [#small[d]:d in [4..6]] eq [126,648,2080];
sid:=0;surv:=[* *];details:=[* *];stats:=[* *];
procedure Test(rows,pattern,~sid,~surv,~details,~trials,~accepts)
 sid+:=1;trials+:=1;ok:=true;
 if exists{r:r in rows|#r[4] ne 2} then ok:=false;else for j in [1,2] do if not &+[q(Z!r[4][j]):r in rows] in H then ok:=false;break;end if;end for;end if;
 if ok then accepts+:=1;Append(~surv,[*sid,[[*r[1],r[3]*]:r in rows]*]);Append(~details,[*sid,pattern,[[r[1],r[2]]:r in rows],[r[4]:r in rows]*]);end if;
end procedure;
for d in [4..10] do trials:=0;accepts:=0;for r in own[d] do Test([*r*],[d],~sid,~surv,~details,~trials,~accepts);end for;Append(~stats,[*[d],trials,accepts*]);print "C242_FIRST_STAGE_PATTERN",[d],trials,accepts;end for;
for pat in [[4,4],[4,5],[4,6],[5,5]] do trials:=0;accepts:=0;for r in own[pat[1]],s in small[pat[2]] do Test([*r,s*],pat,~sid,~surv,~details,~trials,~accepts);end for;Append(~stats,[*pat,trials,accepts*]);print "C242_FIRST_STAGE_PATTERN",pat,trials,accepts;end for;
assert sid eq 7663+n10+25934;
Write(root cat "first_survivors_unified.m",Sprint([*fields,surv*]):Overwrite:=true);Write(root cat "first_survivors_details.m",Sprint([*stats,details*]):Overwrite:=true);print "C242_FIRST_STAGE_COMPLETE",sid,#surv,Cputime(tt);quit;
