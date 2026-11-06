SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(3*10^9);
a:=eval Read("work/deg10/exclusions/315d5/diagonal_data.txt");G4:=AbelianGroup(a[3]);G7:=AbelianGroup(a[5]);Free:=FreeAbelianGroup(20);assert [Eltseq(v):v in Rows(RelationMatrix(G4))] eq a[4];
h4:=hom<Free->G4|[G4!v:v in a[7]]>;H4:=Image(h4);unknown:=[G4!v:v in a[13]];
inside:=[v:v in unknown|v in H4];outside:=[v:v in unknown|v notin H4];print "DECIC315_UNKNOWN_CONTAINMENT",#inside,#outside;
KB:=Matrix(Integers(),[Eltseq(Free!v):v in Generators(Kernel(h4))]);LL:=LatticeWithBasis(LLL(KB));vectors:=[];
for cl in inside do
 v:=Eltseq(cl@@h4);w:=ClosestVector(LL,Vector(Rationals(),v));short:=[v[i]-Integers()!w[i]:i in [1..20]];assert h4(Free!short) eq cl;Append(~vectors,short);
end for;
Write("work/deg10/exclusions/315d5/unknown_cuspidal_lifts.txt",Sprint([*a[1],a[2],vectors,[Eltseq(v):v in outside]*]));print "DECIC315_UNKNOWN_LIFTS_COMPLETE",#vectors;quit;
