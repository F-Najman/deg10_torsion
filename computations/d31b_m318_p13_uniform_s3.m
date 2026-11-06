/* Fresh exact nine-multiple cusp-diagonal inputs; no completion claim. */
SetColumns(0);SetQuitOnError(true);SetNthreads(1);Alarm(1700);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
load "work/deg10/continuation_20261031b/M3N18/inputs/class_group_direct_RR_m3_18_root.m";
Attach("work/deg10/continuation_20261031b/M3N18/inputs/fast_diamond_hash_m2.m");
Attach("work/deg10/continuation_20261029/deps/fast_hecke_11.m");
Attach("work/deg10/continuation_20261031b/M3N18/inputs/fast_full_aut_hash_m3_18_root.m");
for p in [13] do
 X:=MDX11(3,18,GF(p));C:=Curve(X);FF:=FunctionField(C);F0,afmap:=AlgorithmicFunctionField(FF);assert Genus(X) eq 10;
 E:=EllipticCurve([FF!f:f in X`_E]);P:=E![FF!f:f in X`_P];Q:=E![FF!f:f in X`_Q];
 ai:=aInvariants(E);fs:=[ai[3],ai[1]^3-27*ai[3]];
 for k in [1..9] do
  for j in [0..2] do
   A0:=k*Q+j*P;assert A0 ne E!0 and A0[1] ne 0 and A0[2] ne 0;fs cat:=[A0[1],A0[2]];
  end for;
  A0:=k*Q;B0:=A0+P;assert A0 ne B0 and A0 ne -B0;Append(~fs,A0[1]-B0[1]);
 end for;
 cs:=Cusps(X);assert {*Degree(c):c in cs*} eq {*1^^18,3^^6*};groups:=AssociativeArray();
 for c in cs do
  sig:=[Valuation(f,c):f in fs];if IsDefined(groups,sig) then Append(~groups[sig],c);else groups[sig]:=[c];end if;
 end for;
 sigs:=Sort(Setseq(Keys(groups)));assert #sigs eq 24;blocks:=[&+[Divisor(c):c in groups[s]]:s in sigs];ds:=[Degree(b):b in blocks];assert {*d:d in ds*} eq {*1^^18,3^^6*};
 c0:=Support(blocks[Index(ds,1)])[1];
 print "M3N18_9MULT_CUSP_PROFILES",p,#sigs,ds,Cputime(t0);
 G,unused,phi:=ClassGroup(F0:ReductionDivisor:=FunctionFieldDivisor(Divisor(c0)),Proof:=true);
 print "M3N18_9MULT_NATIVE_GROUP_PROVEN",p,Invariants(G),Cputime(t0);
 GG,gm,mg:=MDClassGroupDirectRRCheckedV3(F0,Sprintf("work/deg10/continuation_20261031b/M3N18/paired_uniform_s3_p%o_direct_RR_identities.txt",p));assert GG eq G;
 m2:=map<DivisorGroup(C)->G|D:->mg(FunctionFieldDivisor(D))>;
 gs:=[m2(blocks[i]-ds[i]*Divisor(c0)):i in [1..#blocks]];H:=sub<G|gs>;Qs:=quo<G|H>;
 print "M3N18_9MULT_CUSP_GROUP",p,Invariants(H),Invariants(Qs),Cputime(t0);
 J:=JH(54,[19]);fp:=FrobeniusPolynomial(J,p);assert &*[n:n in Invariants(G)|n ne 0] eq Integers()!Evaluate(fp,1);
 rat:=[c:c in cs|Degree(c) eq 1];ratcl:=[m2(Divisor(c)-Divisor(c0)):c in rat];
 data:=[*p,sigs,ds,Invariants(G),[Eltseq(v):v in Rows(RelationMatrix(G))],[Eltseq(v):v in gs],Invariants(H),[Eltseq(v):v in ratcl],Coefficients(fp)*];
 Write(Sprintf("work/deg10/continuation_20261031b/M3N18/paired_uniform_s3_p%o_group_data.m",p),Sprint(data));

 d:=5;KS:=GF(p^d);PP:=PolynomialRing(KS);qq:=afmap(X`_coordinates[1]);tt:=afmap(X`_coordinates[2]);ypol:=DefiningPolynomial(F0);O:=MaximalOrderFinite(F0);xx:=F0!BaseField(F0).1;yi:=F0.1;ac0:=FunctionFieldDivisor(Divisor(c0));seen:=AssociativeArray();reps:=[* *];vals:=[];portable:=[];images:=[* *];encountered:={};cnt:=0;
 function EvRat(f,b) return Evaluate(Numerator(f),b)/Evaluate(Denominator(f),b);end function;
 function EvFun(f,b,c) return &+[EvRat(Eltseq(f)[i],b)*c^(i-1):i in[1..#Eltseq(f)]];end function;
 function Canon(v) return M3EighteenRootKey(v,KS,p,5);end function;
 function Orbit(q,t)
  EC:=EllipticCurve([q*t-q+t+2,KS!0,q*t^2-q*t+t,KS!0,KS!0]);P0:=EC![0,0];Q0:=EC![-t,t^2];assert Order(P0) eq 3 and Order(Q0) eq 18 and WeilPairing(P0,6*Q0,3) eq KS!RootOfUnity(3,GF(p));ims:=[];
  for par in M3EighteenRootAutParameters() do a,b,c,d1:=Explode(par);Append(~ims,Canon(M3EighteenRootNormalizedPair(EC,a*P0+6*b*Q0,c*P0+d1*Q0)));end for;assert #Seqset(ims) eq 27;return ims;
 end function;
 for beta in KS do
  bc:=[Integers()!a:a in Eltseq(beta,GF(p))];orb:=[bc];b1:=beta;for j in[1..4] do b1:=b1^p;Append(~orb,[Integers()!a:a in Eltseq(b1,GF(p))]);end for;
  if bc ne Sort(orb)[1] then continue;end if;
  if exists{a:a in Coefficients(ypol)|Evaluate(Denominator(a),beta) eq 0} then continue;end if;
  fy:=PP![EvRat(a,beta):a in Coefficients(ypol)];for rr in Roots(fy) do
   gamma:=rr[1];if exists{a:a in Eltseq(qq) cat Eltseq(tt)|Evaluate(Denominator(a),beta) eq 0} then continue;end if;
   qv:=EvFun(qq,beta,gamma);tv:=EvFun(tt,beta,gamma);if qv^p eq qv and tv^p eq tv then continue;end if;
   aa1:=qv*tv-qv+tv+2;aa3:=qv*tv^2-qv*tv+tv;if aa3 eq 0 or aa1^3-27*aa3 eq 0 then continue;end if;
   ky:=Canon([qv,tv]);Include(~encountered,ky);if IsDefined(seen,ky) then continue;end if;
   ims:=Orbit(qv,tv);keep:=(&+[17^j*Sort(ims)[1][j]:j in[1..10]] mod 1000003) mod 4 eq 3;h:=keep select #reps+1 else 0;for key in ims do assert not IsDefined(seen,key);seen[key]:=h;end for;if not keep then continue;end if;
   hb:=MinimalPolynomial(beta,GF(p));if Degree(hb) eq 5 then
    M:=Matrix(GF(p),[Eltseq(beta^j,GF(p)):j in[0..4]]);co:=Eltseq(Solution(M,Vector(GF(p),Eltseq(gamma,GF(p)))));gg:=PolynomialRing(GF(p))!co;I:=ideal<O|Evaluate(hb,xx),yi-Evaluate(gg,xx)>;fac:=Factorization(I);ps:=[Place(r[1]):r in fac|r[2] gt 0 and Degree(Place(r[1])) eq 5];if #ps ne 1 then print "M318_STREAM_IDEAL_DIAGNOSTIC",bc,Eltseq(gamma),hb,co,[(<Degree(r[1]),Degree(Place(r[1])),r[2]>):r in fac];end if;assert #ps eq 1;pl:=ps[1];cu:=CurvePlace(C,pl);
   else
    cu:=MDPlaceCommonZeros(X`_coordinates,[qv,tv]);assert Degree(cu) eq 5;pl:=Support(FunctionFieldDivisor(Divisor(cu)))[1];
   end if;
   assert M3EighteenRootKey([Evaluate(f,cu):f in X`_coordinates],KS,p,5) eq ky;
   cl:=mg(Divisor(pl)-5*ac0);Append(~vals,Eltseq(cl));Append(~portable,[Eltseq(qv),Eltseq(tv)]);Append(~images,ims);Append(~reps,[*bc,Eltseq(gamma)*]);
   if #reps mod 50 eq 0 then Write("work/deg10/continuation_20261031b/M3N18/paired_uniform_s3_partial.m",Sprint([*vals,portable,images,reps,#Keys(seen),#encountered*]):Overwrite:=true);print "M318_STREAM_3_CLASSES",3,#reps,#Keys(seen),#encountered,Cputime(t0);end if;
  end for;
  cnt+:=1;if cnt mod 5000 eq 0 then print "M318_STREAM_BASEPOINTS",cnt,#reps,#Keys(seen),Cputime(t0);end if;
 end for;
 allkeys:=Sort(Setseq(Keys(seen)));assert #allkeys eq 75384;keys:=[v:v in allkeys|seen[v] gt 0];labels:=[seen[v]:v in keys];assert #keys eq 27*#reps;assert #encountered eq #allkeys;
 Write("work/deg10/continuation_20261031b/M3N18/paired_uniform_s3_p13_degree5_data.m",Sprint([*p,5,Coefficients(DefiningPolynomial(KS)),vals,portable,[1..#reps],keys,labels,images*]):Overwrite:=true);
 print "D31B_M318_UNIFORM_P13_3_COMPLETE",3,#keys,#reps,Cputime(t0);
end for;
quit;
