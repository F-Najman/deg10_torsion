SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(3*10^9);
a:=eval Read("work/deg10/exclusions/315d5/diagonal_data.txt");G4:=AbelianGroup(a[3]);G7:=AbelianGroup(a[5]);Free:=FreeAbelianGroup(20);h4:=hom<Free->G4|[G4!v:v in a[7]]>;h7:=hom<Free->G7|[G7!v:v in a[8]]>;H4:=Image(h4);H7:=Image(h7);
known:={G4|G4!v:v in a[12]};surv:={G4|};pairs:=[];
for p in a[11] do
 cl4:=G4!a[9][p[1]][3];cl7:=G7!a[10][p[2]][3];
 if cl4 notin H4 then continue;end if;
 if h7(Free!(cl4@@h4)) eq cl7 then Include(~surv,cl4);Append(~pairs,p);end if;
end for;
unknown:=surv diff known;print "DECIC315_EXACT_DIAGONAL",#pairs,#surv,#unknown;
Write("work/deg10/exclusions/315d5/exact_pairs.txt",Sprint([*pairs,[Eltseq(v):v in surv],[Eltseq(v):v in unknown]*]));quit;
