SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(5*10^9);Alarm(800);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
Attach("work/deg10/continuation_20261029/deps/fast_hecke.m");
Attach("work/deg10/continuation_20261029/deps/fast_hecke_11.m");
load "work/deg10/continuation_20261031b/M3N18/fast_class_group_v4.m";
Attach("work/deg10/continuation_20261031b/M3N18/inputs/fast_full_aut_hash_m3_18_root.m");
X:=MDX11(3,18,GF(7));C:=Curve(X);FF:=FunctionField(C);E:=EllipticCurve(X`_E);P:=E!X`_P;Q:=E!X`_Q;
ai:=aInvariants(E);fs:=[ai[3],ai[1]^3-27*ai[3]];
for k in [1..9] do
 for j in [0..2] do A:=k*Q+j*P;assert A ne E!0 and A[1] ne 0 and A[2] ne 0;fs cat:=[A[1],A[2]];end for;
 A:=k*Q;B:=A+P;assert A ne B and A ne -B;Append(~fs,A[1]-B[1]);
end for;
cs:=Cusps(X);groups:=AssociativeArray();
for c in cs do
 sig:=[Valuation(f,c):f in fs];if IsDefined(groups,sig) then Append(~groups[sig],c);else groups[sig]:=[c];end if;
end for;
sigs:=Sort(Setseq(Keys(groups)));blocks:=[&+[Divisor(c):c in groups[sig]]:sig in sigs];ds:=[Degree(b):b in blocks];assert #blocks eq 24;
c0:=Support(blocks[Index(ds,1)])[1];AF:=AlgorithmicFunctionField(FF);
G,_,native:=ClassGroup(AF:ReductionDivisor:=FunctionFieldDivisor(Divisor(c0)),Proof:=true);
GG,_,ph:=MDClassGroupDirectRRFastV4(AF,"work/deg10/continuation_20261031b/M3N18/inert5_framed_RR_identities.log");assert GG eq G;
m2:=map<DivisorGroup(C)->G|D:->ph(FunctionFieldDivisor(D))>;
cc:=[m2(blocks[i]-ds[i]*Divisor(c0)):i in [1..24]];Hc:=sub<G|cc>;L:=Hc;
order:=&*[e:e in Invariants(G)|e ne 0];extra:=[G|];picked:=[* *];
for d in [2,3,4] do
 if #L eq order then break;end if;xs:=NoncuspidalPlaces(X,d);
 for i in [1..#xs] do
  cl:=m2(Divisor(xs[i])-d*Divisor(c0));
  if cl notin L then Append(~extra,cl);Append(~picked,[*d,i*]);L:=sub<G|Setseq(Generators(L)) cat [cl]>;
   print "DECIC_HECKE_TORSION_ADD",3,18,d,i,#L div #Hc,Cputime(tt);end if;
  if #L eq order then break;end if;
 end for;
end for;
portable:=[* *];for pair in picked do dd:=pair[1];pl:=NoncuspidalPlaces(X,dd)[pair[2]];KK:=GF(7^dd);qt:=M3EighteenRootToCommon([Evaluate(v,pl):v in X`_coordinates],KK);Append(~portable,[*dd,Coefficients(DefiningPolynomial(KK)),[Eltseq(v):v in qt]*]);end for;Write("work/deg10/continuation_20261031b/M3N18/inert5_generator_points.m",Sprint(portable):Overwrite:=true);
old:=eval Read("work/deg10/continuation_20261031b/M3N18/split_hecke_framed_arithmetic.txt");assert [Eltseq(c):c in cc]eq old[5] and picked eq old[6] and [Eltseq(c):c in extra]eq old[7];assert #L eq order;Free:=FreeAbelianGroup(24+#extra);h:=hom<Free->G|cc cat extra>;GG:=Image(h);assert #GG eq order;K:=GG;rows:=[* *];
print "DECIC_HECKE_TORSION_SETUP",3,18,#Hc,Invariants(quo<GG|sub<GG|[GG!c:c in cc]>>),picked,Cputime(tt);
for q in [5] do
 imgs:=[G!0:i in [1..24]];
 for pair in picked do
  d:=pair[1];x:=NoncuspidalPlaces(X,d)[pair[2]];
  T:=HeckeOperatorQuotientFast(X,5,x);A:=HeckeOperatorQuotientFast(X,5,T)-Divisor(x)-10*Divisor(DiamondOperator(X,5,x))-25*Divisor(DiamondOperator(X,25,x));
  assert Degree(A) eq 0;Append(~imgs,m2(A));
 end for;
 for v in Generators(Kernel(h)) do assert &+[Eltseq(Free!v)[i]*imgs[i]:i in [1..#imgs]] eq G!0;end for;
 out:=[&+[Eltseq((G!GG.j)@@h)[i]*imgs[i]:i in [1..#imgs]]:j in [1..Ngens(GG)]];
 aq:=hom<GG->G|out>;K meet:=Kernel(aq);assert #K mod #Hc eq 0;
 Append(~rows,[*q,[Eltseq(g):g in imgs],[Eltseq(g):g in out],#K div #Hc*]);
 print "DECIC_HECKE_TORSION_KERNEL",3,18,q,#K div #Hc,Cputime(tt);
 if #K eq #Hc then break;end if;
end for;
Write("work/deg10/continuation_20261031b/M3N18/inert5_framed_arithmetic.m",Sprint([*sigs,ds,Invariants(G),[Eltseq(r):r in Rows(RelationMatrix(G))],[Eltseq(c):c in cc],picked,[Eltseq(c):c in extra],Invariants(GG),rows,#K div #Hc*]));
print "D31B_M318_INERT5_ENDOMORPHISMS_COMPLETE",3,18,#K div #Hc,Cputime(tt);quit;
