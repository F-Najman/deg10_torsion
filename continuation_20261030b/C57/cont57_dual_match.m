SetColumns(0);SetNthreads(1);SetQuitOnError(true);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
load "work/deg10/continuation_20261030b/C57/census10_records.m";assert N eq 57 and FieldDegree eq 10;
F<a>:=ext<GF(2)|PolynomialRing(GF(2))!FieldPolynomial>;
ev:=func<v|Evaluate(PolynomialRing(GF(2))!v,a)>;
rs:=[r:r in Records];ivs:=SequenceToMultiset([r[5]:r in rs]);print "DECIC57_POINT_GROUPS",ivs;
data:=eval Read("work/deg10/continuation_20261030b/C57/noncuspidal_filter_v2_data.txt");kept:=data[2];residual_jps:={r[4]:r in kept};assert #residual_jps eq 3;jps:={r[6]:r in rs};assert #jps eq 6;
Code:=func<v|[Integers()!c:c in Eltseq(v,GF(2))]>;
function PKey(P)
 return Code(P[1]) cat Code(P[2]);
end function;
load "work/deg10/continuation_20261030b/C57/selected_levels.m";previous:=SelectedLevels;selected:=[* *];all_rows:=[* *];index:=0;missing:=[* *];
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
  RR<x>:=PolynomialRing(F);kp:=&*[x-(u*Q)[1]:u in [1..28]];assert Degree(kp) eq 28;EQ,iso:=IsogenyFromKernel(E,kp);assert Degree(iso) eq 57;dualjp:=[Integers()!v:v in Eltseq(MinimalPolynomial(jInvariant(EQ),GF(2)))];needed:=dualjp in residual_jps;
  old:=[v[1]:v in previous|v[2] eq jp and v[6] eq key];assert #old le 1;
  rec:=[*index,jp,iv,[Code(v):v in aInvariants(E)],[Code(Q[1]),Code(Q[2])],key*];Append(~all_rows,[*rec,dualjp,needed,old*]);
  if needed then Append(~selected,rec);if #old eq 0 then Append(~missing,rec);end if;end if;
  printf "DECIC57_PHYSICAL_DUAL_J_MATCH source=%o dual=%o needed=%o prior_test=%o\n",jp,dualjp,needed,old;
 end for;
 assert #allP eq 36*num;printf "DECIC57_SOURCE_LEVEL_ORBITS jpoly=%o group=%o primitive_points=%o diamond_orbits=%o\n",jp,iv,#allP,num;
end for;
Write("work/deg10/continuation_20261030b/C57/matched_levels.m",Sprintf("FieldPolynomial:=%o;MatchedLevels:=%o;MissingLevels:=%o;AllSourceMatches:=%o;",FieldPolynomial,selected,missing,all_rows):Overwrite:=true);
printf "DECIC57_PHYSICAL_DUAL_MATCH_COMPLETE cases=%o missing=%o time=%o\n",#selected,#missing,Cputime(tt);quit;
