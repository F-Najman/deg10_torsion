SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(4*10^9);tt:=Cputime();
root:="work/deg10/C236_resolution_20261030/enumeration/";
ad:=eval Read(root cat "target18_aut_actions.m");Z:=FreeAbelianGroup(#ad[2][1]);Q,q:=quo<Z|[Z!r:r in ad[2]]>;assert Invariants(Q) eq [2,6];
function RawClass(v) assert #v eq 4 and v[4] eq 0;return &+[v[j]*(Q!ad[3][j]):j in [1..3]];end function;
own:=AssociativeArray();small:=AssociativeArray();fields:=[* *];
for d in [3..7] do
 a:=eval Read(root cat Sprintf("small_qclasses%o.m",d));Append(~fields,[*d,a[1]*]);own[d]:=a[2];small[d]:=a[3];
end for;
for d in [8,9] do
 a:=eval Read(root cat Sprintf("filter18_D%o_D%o.m",d,d));Append(~fields,[*d,a[1]*]);own[d]:=[*[*r[1],r[2],r[3],[Eltseq(RawClass(v)):v in r[4]]*]:r in a[2]*];
end for;
a:=eval Read(root cat "filter18_D10_D10s0.m");field10:=a[1];dict:=AssociativeArray();
for s in [0..7] do
 a:=eval Read(root cat Sprintf("filter18_D10_D10s%o.m",s));assert a[1] eq field10;
 for r in a[2] do assert not IsDefined(dict,r[2]);dict[r[2]]:=[*r[1],r[2],r[3],[Eltseq(RawClass(v)):v in r[4]]*];end for;
end for;
assert Sort(Setseq(Keys(dict))) eq [1..40874];own[10]:=[*dict[i]:i in [1..40874]*];Append(~fields,[*10,field10*]);
assert [#own[d]:d in [3..10]] eq [2,5,30,126,434,2115,9088,40874];
for d in [3..10] do assert [r[2]:r in own[d]] eq [1..#own[d]];end for;
sid:=0;surv:=[* *];details:=[* *];stats:=[* *];
procedure Test(rows,pattern,~sid,~surv,~details,~trials,~accepts)
 sid+:=1;trials+:=1;ok:=true;
 if exists{r:r in rows|#r[4] ne 2} then ok:=false;else
  for j in [1,2] do if &+[Q!r[4][j]:r in rows] ne Q!0 then ok:=false;break;end if;end for;
 end if;
 if ok then accepts+:=1;Append(~surv,[*sid,[[*r[1],r[3]*]:r in rows]*]);Append(~details,[*sid,pattern,[[r[1],r[2]]:r in rows],[r[4]:r in rows]*]);end if;
end procedure;
for d in [3..10] do
 trials:=0;accepts:=0;for r in own[d] do Test([*r*],[d],~sid,~surv,~details,~trials,~accepts);end for;
 Append(~stats,[*[d],trials,accepts*]);print "C236_FIRST_STAGE_PATTERN",[d],trials,accepts;
end for;
for pat in [[3,3],[3,4],[3,5],[3,6],[3,7],[4,4],[4,5],[4,6],[5,5]] do
 trials:=0;accepts:=0;for r in own[pat[1]],s in small[pat[2]] do Test([*r,s*],pat,~sid,~surv,~details,~trials,~accepts);end for;
 Append(~stats,[*pat,trials,accepts*]);print "C236_FIRST_STAGE_PATTERN",pat,trials,accepts;
end for;
trials:=0;accepts:=0;for r in own[3] do for j in [1..#small[3]],k in [j..#small[3]] do Test([*r,small[3][j],small[3][k]*],[3,3,3],~sid,~surv,~details,~trials,~accepts);end for;end for;
Append(~stats,[*[3,3,3],trials,accepts*]);print "C236_FIRST_STAGE_PATTERN",[3,3,3],trials,accepts;
trials:=0;accepts:=0;for r in own[3],s in small[3],t in small[4] do Test([*r,s,t*],[3,3,4],~sid,~surv,~details,~trials,~accepts);end for;
Append(~stats,[*[3,3,4],trials,accepts*]);print "C236_FIRST_STAGE_PATTERN",[3,3,4],trials,accepts;
assert sid eq 135282;
Write(root cat "first_survivors_unified.m",Sprint([*fields,surv*]):Overwrite:=true);
Write(root cat "first_survivors_details.m",Sprint([*stats,details*]):Overwrite:=true);
print "C236_FIRST_STAGE_COMPLETE",sid,#surv,Cputime(tt);quit;
