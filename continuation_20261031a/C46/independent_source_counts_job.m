SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(4*10^9);t0:=Cputime();
AttachSpec("work/degree9/shared/sutherland-magma-current/magma.spec");
for tup in [<46,3>] do
 n,p:=Explode(tup);GG:=GL(2,Integers(n));
 HH:=sub<GG|[GG![1,0,1,1]] cat [GG![a,0,0,1]:a in [1..n]|GCD(a,n) eq 1 ]>;
 zs:=GL2PointCounts(HH,p,10);
 ys:=[zs[d]-GL2RationalCuspCount(HH,p^d):d in [1..10]];
 assert &+GL2jCounts(HH,p) eq ys[1];closed:=[];
 for d in [1..10] do
  z:=ys[d]-&+([0] cat [e*closed[e]:e in [1..d-1]|d mod e eq 0]);
  assert z mod d eq 0;Append(~closed,z div d);
 end for;
 assert closed eq [0,0,0,55,44,55,264,781,1958,6292];
 Write("work/deg10/continuation_20261031a/C46/independent_source_counts.m",Sprintf("ClosedCounts:=%o;",closed));
 printf "D31A_C46_INDEPENDENT_COUNT m=1 n=%o p=%o genus=%o rational_point_counts=%o closed_place_counts=%o time=%o\n",n,p,GL2Genus(HH),ys,closed,Cputime(t0);
end for;
print "D31A_C46_INDEPENDENT_COUNTS_DONE complete=true",Cputime(t0);
quit;
