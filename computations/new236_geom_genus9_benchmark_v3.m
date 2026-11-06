a:=eval Read("work/deg10/C236_resolution_20261030/geometry/quotient9_proven_picard.m");a[13]:=a[10];Write("work/deg10/C236_resolution_20261030/geometry/quotient9_proven_picard_v2.m",Sprint(a):Overwrite:=true);
src:=eval Read("work/deg10/C236_resolution_20261030/enumeration/source10_s0.m");fields:=[[*10,src[2]*]];rr:=[r:r in src[9]|r[1] eq 10];cases:=[*[*i,[*[*10,rr[i][2]*]*]*]:i in [1..5]*];Write("work/deg10/C236_resolution_20261030/geometry/genus9_benchmark_v3_input.m",Sprint([*fields,cases*]):Overwrite:=true);
INPUTPATH:="work/deg10/C236_resolution_20261030/geometry/genus9_benchmark_v3_input.m";OUTPUTPREFIX:="work/deg10/C236_resolution_20261030/geometry/genus9_benchmark_v3";MAXCASES:=5;
load "work/deg10/C236_resolution_20261030/geometry/genus9_stage2_common_v2.m";
