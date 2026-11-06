a:=eval Read("work/deg10/C236_resolution_20261030/geometry/quotient9_proven_picard.m");a[13]:=a[10];Write("work/deg10/C236_resolution_20261030/geometry/quotient9_proven_picard_v2.m",Sprint(a):Overwrite:=true);
src:=eval Read("work/deg10/C236_resolution_20261030/enumeration/source10_s0.m");fields:=[[*10,src[2]*]];cases:=[*[*i,[*[*10,src[9][i][2]*]*]*]:i in [1..5]*];Write("work/deg10/C236_resolution_20261030/geometry/genus9_benchmark_v2_input.m",Sprint([*fields,cases*]):Overwrite:=true);
INPUTPATH:="work/deg10/C236_resolution_20261030/geometry/genus9_benchmark_v2_input.m";OUTPUTPREFIX:="work/deg10/C236_resolution_20261030/geometry/genus9_benchmark_v2";MAXCASES:=5;
load "work/deg10/C236_resolution_20261030/geometry/genus9_stage2_common_v2.m";
