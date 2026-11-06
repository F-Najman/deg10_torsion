from pathlib import Path
H=Path(__file__).resolve().parents[1]
for m,n in [(3,18),(6,12)]:
 stem=f'cont_hecke_torsion_{m}_{n}'
 s='''SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(5*10^9);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
Attach("work/deg10/continuation_20261029/deps/fast_hecke.m");
Attach("work/deg10/continuation_20261029/deps/fast_hecke_11.m");
load "work/deg10/continuation_20261029/deps/class_group_direct_RR_m3_18_root.m";
'''
 if m==3:
  s+='''X:=MDX11(3,18,GF(7));C:=Curve(X);FF:=FunctionField(C);E:=EllipticCurve(X`_E);P:=E!X`_P;Q:=E!X`_Q;
ai:=aInvariants(E);fs:=[ai[3],ai[1]^3-27*ai[3]];
for k in [1..9] do
 for j in [0..2] do A:=k*Q+j*P;assert A ne E!0 and A[1] ne 0 and A[2] ne 0;fs cat:=[A[1],A[2]];end for;
 A:=k*Q;B:=A+P;assert A ne B and A ne -B;Append(~fs,A[1]-B[1]);
end for;
'''
 else:
  s+='''Attach("work/deg10/continuation_20261029/deps/model_6_12.m");
X:=MDX612(GF(7));C:=Curve(X);FF:=FunctionField(C);E:=EllipticCurve(X`_E);Q:=E![0,0];px:=FF!X`_P[1];ai:=aInvariants(E);R<Y>:=PolynomialRing(FF);
ys:=Roots(Y^2+(ai[1]*px+ai[3])*Y-px^3-ai[2]*px^2-ai[4]*px-ai[5]);assert #ys eq 2;P:=E![px,ys[1][1]];
z:=PrimitiveElement(GF(7));if WeilPairing(P,2*Q,6) ne z then P:=-P;end if;assert WeilPairing(P,2*Q,6) eq z;
fs:=[-ai[2],1-ai[1],Discriminant(E)];anchor:=P[1];
for i in [0..5] do for j in [0..11] do
 A:=i*P+j*Q;if A eq E!0 then continue;end if;
 if A ne Q and A ne -Q then assert A[1] ne 0;Append(~fs,A[1]);end if;
 if A ne P and A ne -P then assert A[1] ne anchor;Append(~fs,A[1]-anchor);end if;
end for;end for;
'''
 s+=f'''cs:=Cusps(X);groups:=AssociativeArray();
for c in cs do
 sig:=[Valuation(f,c):f in fs];if IsDefined(groups,sig) then Append(~groups[sig],c);else groups[sig]:=[c];end if;
end for;
sigs:=Sort(Setseq(Keys(groups)));blocks:=[&+[Divisor(c):c in groups[sig]]:sig in sigs];ds:=[Degree(b):b in blocks];assert #blocks eq 24;
c0:=Support(blocks[Index(ds,1)])[1];AF:=AlgorithmicFunctionField(FF);
G,_,native:=ClassGroup(AF:ReductionDivisor:=FunctionFieldDivisor(Divisor(c0)),Proof:=true);
GG,_,ph:=MDClassGroupDirectRRCheckedV3(AF,"work/deg10/continuation_20261029/{stem}_RR_identities.log");assert GG eq G;
m2:=map<DivisorGroup(C)->G|D:->ph(FunctionFieldDivisor(D))>;
cc:=[m2(blocks[i]-ds[i]*Divisor(c0)):i in [1..24]];Hc:=sub<G|cc>;L:=Hc;
order:=&*[e:e in Invariants(G)|e ne 0];extra:=[G|];picked:=[* *];
for d in [2,3,4] do
 if #L eq order then break;end if;xs:=NoncuspidalPlaces(X,d);
 for i in [1..#xs] do
  cl:=m2(Divisor(xs[i])-d*Divisor(c0));
  if cl notin L then Append(~extra,cl);Append(~picked,[*d,i*]);L:=sub<G|Generators(L) cat [cl]>;
   print "DECIC_HECKE_TORSION_ADD",{m},{n},d,i,#L div #Hc,Cputime(tt);end if;
  if #L eq order then break;end if;
 end for;
end for;
assert #L eq order;Free:=FreeAbelianGroup(24+#extra);h:=hom<Free->G|cc cat extra>;GG:=Image(h);assert #GG eq order;K:=GG;rows:=[* *];
print "DECIC_HECKE_TORSION_SETUP",{m},{n},#Hc,Invariants(quo<GG|sub<GG|[GG!c:c in cc]>>),picked,Cputime(tt);
for q in [13,19,31,37] do
 imgs:=[G!0:i in [1..24]];
 for pair in picked do
  d:=pair[1];x:=NoncuspidalPlaces(X,d)[pair[2]];
'''
 if m==3:s+='  A:=HeckeOperatorFast(X,q,x)-Divisor(x)-q*Divisor(DiamondOperator(X,q,x));\n'
 else:s+='  A:=HeckeOperator612(X,q,x)-Divisor(x)-q*Divisor(Diamond612(X,q,x));\n'
 s+=f'''  assert Degree(A) eq 0;Append(~imgs,m2(A));
 end for;
 for v in Generators(Kernel(h)) do assert &+[Eltseq(Free!v)[i]*imgs[i]:i in [1..#imgs]] eq G!0;end for;
 out:=[&+[Eltseq((G!GG.j)@@h)[i]*imgs[i]:i in [1..#imgs]]:j in [1..Ngens(GG)]];
 aq:=hom<GG->G|out>;K meet:=Kernel(aq);assert #K mod #Hc eq 0;
 Append(~rows,[*q,[Eltseq(g):g in imgs],[Eltseq(g):g in out],#K div #Hc*]);
 print "DECIC_HECKE_TORSION_KERNEL",{m},{n},q,#K div #Hc,Cputime(tt);
 if #K eq #Hc then break;end if;
end for;
Write("work/deg10/continuation_20261029/{stem}_arithmetic.txt",Sprint([*sigs,ds,Invariants(G),[Eltseq(r):r in Rows(RelationMatrix(G))],[Eltseq(c):c in cc],picked,[Eltseq(c):c in extra],Invariants(GG),rows,#K div #Hc*]));
print "DECIC_HECKE_TORSION_SATURATION_COMPLETE",{m},{n},#K div #Hc,Cputime(tt);quit;
'''
 (H/'computations'/(stem+'.m')).write_text(s)
