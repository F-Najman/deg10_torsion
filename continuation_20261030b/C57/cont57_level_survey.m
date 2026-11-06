SetColumns(0);SetNthreads(1);SetQuitOnError(true);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
load "work/deg10/continuation_20261030b/C57/census10_records.m";assert N eq 57 and FieldDegree eq 10;
F<a>:=ext<GF(2)|PolynomialRing(GF(2))!FieldPolynomial>;
ev:=func<v|Evaluate(PolynomialRing(GF(2))!v,a)>;
rs:=[r:r in Records];ivs:=SequenceToMultiset([r[5]:r in rs]);print "DECIC57_POINT_GROUPS",ivs;
data:=eval Read("work/deg10/continuation_20261030b/C57/noncuspidal_filter_v2_data.txt");kept:=data[2];jps:={r[4]:r in kept};assert #jps eq 3;
Code:=func<v|[Integers()!c:c in Eltseq(v,GF(2))]>;
function PKey(P)
 return Code(P[1]) cat Code(P[2]);
end function;
selected:=[* *];index:=0;
for jp in Sort(Setseq(jps)) do
 rows:=[r:r in rs|r[6] eq jp];assert #rows eq 10;r:=rows[1];E:=EllipticCurve([ev(v):v in r[2]]);A,mp:=AbelianGroup(E);iv:=Invariants(A);allP:=[];
 if #iv eq 1 then
  for k in [0..GCD(57,iv[1])-1] do Q:=mp((k*(iv[1] div GCD(57,iv[1])))*A.1);if Order(Q) eq 57 then Append(~allP,Q);end if;end for;
 else
  assert #iv eq 2;
  for k in [0..GCD(57,iv[1])-1],l in [0..GCD(57,iv[2])-1] do Q:=mp(k*(iv[1] div GCD(57,iv[1]))*A.1+l*(iv[2] div GCD(57,iv[2]))*A.2);if Order(Q) eq 57 then Append(~allP,Q);end if;end for;
 end if;
 seen:={};num:=0;
 for Q in allP do
  key:=Min([PKey(u*Q):u in [1..56]|GCD(u,57) eq 1]);if key in seen then continue;end if;Include(~seen,key);num+:=1;index+:=1;
  rec:=[*index,jp,iv,[Code(v):v in aInvariants(E)],[Code(Q[1]),Code(Q[2])],key*];Append(~selected,rec);
 end for;
 assert #allP eq 36*num;printf "DECIC57_SOURCE_LEVEL_ORBITS jpoly=%o group=%o primitive_points=%o diamond_orbits=%o\n",jp,iv,#allP,num;
end for;
Write("work/deg10/continuation_20261030b/C57/selected_levels.m",Sprintf("FieldPolynomial:=%o;SelectedLevels:=%o;",FieldPolynomial,selected):Overwrite:=true);
printf "DECIC57_LEVEL_SURVEY_COMPLETE cases=%o time=%o\n",#selected,Cputime(tt);quit;
