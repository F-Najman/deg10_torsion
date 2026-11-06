SetColumns(0);SetQuitOnError(true);SetNthreads(1);tt:=Cputime();
a:=eval Read("work/deg10/continuation_20261029/diagonal612_D5_source_checkpoint_v5.m");
FB:=FreeAbelianGroup(#a[5][2][1]);G,pi:=quo<FB|[FB!r:r in a[4][2]]>;Free:=FreeAbelianGroup(24);h:=hom<Free->G|[pi(FB!v):v in a[5][2]]>;assert #Image(h) eq 5971968;
KB:=Matrix(Integers(),[Eltseq(Free!w):w in Generators(Kernel(h))]);LL:=LatticeWithBasis(LLL(KB));classes:=AssociativeArray();lifts:=[* *];
for i in [1..#a[8]] do
 row:=a[8][i];if not row[5] then continue;end if;key:=row[6];
 if IsDefined(classes,key) then Append(~lifts[classes[key]][1],i);continue;end if;
 cl:=pi(FB!key);vv:=Eltseq(cl@@h);w:=ClosestVector(LL,Vector(Rationals(),vv));v:=[vv[j]-Integers()!w[j]:j in [1..24]];assert h(Free!v) eq cl;
 Append(~lifts,[*[i],v,key*]);classes[key]:=#lifts;
end for;
assert #lifts eq 3360;Write("work/deg10/continuation_20261029/M6N12_D5_all_source_lifts.m",Sprint([*a[1],a[2],lifts*]));print "DECIC612_ALL_SOURCE_CUSP_LIFTS_COMPLETE",#lifts,Cputime(tt);quit;
