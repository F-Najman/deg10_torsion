SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(8*10^9);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
Attach("work/deg10/continuation_20261029/deps/fast_hecke.m");
Attach("work/deg10/continuation_20261029/deps/fast_hecke_11.m");
Attach("work/deg10/continuation_20261029/deps/fast_diamond_hash_m2.m");
Attach("work/deg10/continuation_20261029/deps/model_6_12.m");
Attach("work/deg10/continuation_20261029/deps/full_aut_6_12.m");
Attach("work/deg10/continuation_20261029/deps/full_aut_6_12_compact.m");
load "work/deg10/continuation_20261029/deps/class_group_direct_RR_m3_18_root.m";
function CuspBlocks(X)
C:=Curve(X);FF:=FunctionField(C);E:=EllipticCurve(X`_E);Q:=E![0,0];px:=FF!X`_P[1];ai:=aInvariants(E);R<Y>:=PolynomialRing(FF);
ys:=Roots(Y^2+(ai[1]*px+ai[3])*Y-px^3-ai[2]*px^2-ai[4]*px-ai[5]);assert #ys eq 2;P:=E![px,ys[1][1]];
z:=PrimitiveElement(BaseRing(X))^((#BaseRing(X)-1) div 6);if WeilPairing(P,2*Q,6) ne z then P:=-P;end if;assert WeilPairing(P,2*Q,6) eq z;
fs:=[-ai[2],1-ai[1],Discriminant(E)];anchor:=P[1];
for i in [0..5] do for j in [0..11] do
 A:=i*P+j*Q;if A eq E!0 then continue;end if;
 if A ne Q and A ne -Q then assert A[1] ne 0;Append(~fs,A[1]);end if;
 if A ne P and A ne -P then assert A[1] ne anchor;Append(~fs,A[1]-anchor);end if;
end for;end for;
cs:=Cusps(X);groups:=AssociativeArray();
for c in cs do
 sig:=[Valuation(f,c):f in fs];if IsDefined(groups,sig) then Append(~groups[sig],c);else groups[sig]:=[c];end if;
end for;
sigs:=Sort(Setseq(Keys(groups)));blocks:=[&+[Divisor(c):c in groups[sig]]:sig in sigs];ds:=[Degree(b):b in blocks];assert #blocks eq 24;
return blocks,sigs,ds;
end function;
function WeightedVectors(ds,d)
 if #ds eq 0 then return d eq 0 select [[]] else [];end if;
 out:=[];for k in [0..d div ds[1]] do for v in WeightedVectors(ds[2..#ds],d-k*ds[1]) do Append(~out,[k] cat v);end for;end for;return out;
end function;
Free:=FreeAbelianGroup(24);Xs:=[* *];Gs:=[* *];hs:=[* *];ms:=[* *];bases:=[* *];cgs:=[* *];sigs:=[];ds:=[];
for p in [7,13] do
 X:=MDX612(GF(p));C:=Curve(X);assert Genus(X) eq 9;blocks,ss,dd:=CuspBlocks(X);
 if p eq 7 then sigs:=ss;ds:=dd;else assert ss eq sigs and dd eq ds;end if;
 c0:=Support(blocks[Index(ds,1)])[1];AF:=AlgorithmicFunctionField(FunctionField(C));
 G,_,native:=ClassGroup(AF:ReductionDivisor:=FunctionFieldDivisor(Divisor(c0)),Proof:=true);
 G,_,ph:=MDClassGroupDirectRRCheckedV3(AF,Sprintf("work/deg10/continuation_20261029/diagonal612_D5_v4_RR%o_identities.log",p));
 m2:=map<DivisorGroup(C)->G|D:->ph(FunctionFieldDivisor(D))>;
 cc:=[m2(blocks[i]-ds[i]*Divisor(c0)):i in [1..24]];h:=hom<Free->G|cc>;assert #Image(h) eq 5971968;
 Append(~Xs,X);Append(~Gs,G);Append(~hs,h);Append(~ms,m2);Append(~bases,c0);Append(~cgs,cc);
 print "DECIC612_D5_GROUP",p,Invariants(G),Invariants(Image(h)),Cputime(tt);
end for;
h7:=hs[1];h13:=hs[2];G7:=Gs[1];G13:=Gs[2];
for v in Generators(Kernel(h7)) do assert h13(Free!v) eq G13!0;end for;
for v in Generators(Kernel(h13)) do assert h7(Free!v) eq G7!0;end for;
X:=Xs[1];C:=Curve(X);m2:=ms[1];c0:=bases[1];L:=Image(h7);extra:=[G7|];picked:=[];xs:=NoncuspidalPlaces(X,3);assert #xs eq 96;
for i in [1..#xs] do
 cl:=m2(Divisor(xs[i])-3*Divisor(c0));
 if cl notin L then Append(~extra,cl);Append(~picked,i);L:=sub<G7|Setseq(Generators(L)) cat [cl]>;end if;
 if #L eq 95551488 then break;end if;
end for;
assert #L eq 95551488;F:=FreeAbelianGroup(24+#extra);h:=hom<F->G7|cgs[1] cat extra>;GG:=Image(h);
imgs:=[G7!0:i in [1..24]];
for i in picked do
 x:=xs[i];A:=HeckeOperator612(X,13,x)-Divisor(x)-13*Divisor(Diamond612(X,13,x));Append(~imgs,m2(A));
end for;
for v in Generators(Kernel(h)) do assert &+[Eltseq(F!v)[i]*imgs[i]:i in [1..#imgs]] eq G7!0;end for;
out:=[&+[Eltseq((G7!GG.j)@@h)[i]*imgs[i]:i in [1..#imgs]]:j in [1..Ngens(GG)]];
aq:=hom<GG->G7|out>;K:=Kernel(aq);assert Invariants(quo<K|sub<K|[K!c:c in cgs[1]]>>) eq [2,2];
print "DECIC612_D5_HECKE_KERNEL",Invariants(K),Invariants(quo<K|sub<K|[K!c:c in cgs[1]]>>),Cputime(tt);
old:=eval Read("work/deg10/continuation_20261029/diagonal612_D5_source_checkpoint.m");assert sigs eq old[1] and ds eq old[2];
OB:=FreeAbelianGroup(#old[5][1][1]);OG,op:=quo<OB|[OB!r:r in old[4][1]]>;
occ:=[op(OB!v):v in old[5][1]];oldextra:=[OG|];newextra:=[G7|];
assert picked eq [1,2];
for i in picked do
 row:=[r:r in old[8]|r[1] eq 3 and r[2] eq i][1];
 Append(~oldextra,op(OB!row[4])-&+[row[3][j]*occ[j]:j in [1..24]]);
 Append(~newextra,m2(Divisor(xs[i])-3*Divisor(c0)));
end for;
OF:=FreeAbelianGroup(26);oh:=hom<OF->OG|occ cat oldextra>;nh:=hom<OF->G7|cgs[1] cat newextra>;
assert #Image(oh) eq 95551488 and #Image(nh) eq 95551488;
for v in Generators(Kernel(oh)) do assert nh(OF!v) eq G7!0;end for;
for v in Generators(Kernel(nh)) do assert oh(OF!v) eq OG!0;end for;
source:=[* *];predictions:=AssociativeArray();cache:=AssociativeArray();sourcecoverage:=old[9];skipped:=0;
for row in old[8] do
 d:=row[1];i:=row[2];v:=row[3];ocl:=op(OB!row[4]);cl:=nh(OF!(ocl@@oh));pass:=cl in K;assert pass eq row[5];
 if not pass then Append(~source,[*d,i,v,Eltseq(cl),false,[]*]);skipped+:=1;continue;end if;
 tw:=2*cl;assert tw in Image(h7);key:=Eltseq(tw);
 if not IsDefined(cache,key) then cache[key]:=h13(Free!(tw@@h7));end if;
 pred:=cache[key];Append(~source,[*d,i,v,Eltseq(cl),true,Eltseq(pred)*]);ky:=Eltseq(pred);
 if IsDefined(predictions,ky) then Append(~predictions[ky],#source);else predictions[ky]:=[#source];end if;
end for;
assert #source eq 24576;print "DECIC612_D5_SOURCE_RESUME_COMPLETE",#source,skipped,#Keys(predictions),Cputime(tt);
prefix:=[*sigs,ds,[Invariants(G):G in Gs],[[Eltseq(r):r in Rows(RelationMatrix(G))]:G in Gs],[[Eltseq(c):c in cg]:cg in cgs],[Eltseq(G7!GG.j):j in [1..Ngens(GG)]],[Eltseq(g):g in out],source,sourcecoverage*];
Write("work/deg10/continuation_20261029/diagonal612_D5_source_checkpoint_v4.m",Sprint(prefix));
X:=Xs[2];C:=Curve(X);m2:=ms[2];c0:=bases[2];rawcs:=Cusps(X);rawds:=[Degree(c):c in rawcs];rc:=[m2(Divisor(c)-Degree(c)*Divisor(c0)):c in rawcs];
targets:=[* *];matches:=[* *];coverage:=[* *];counts:=[];allxs:=[* *];allcl:=[* *];reps:=[* *];repc:=[* *];
procedure Check(cl,descr,~targets,~matches)
 Append(~targets,[*descr,Eltseq(cl)*]);ky:=Eltseq(2*cl);
 if IsDefined(predictions,ky) then Append(~matches,[*#targets,predictions[ky]*]);end if;
end procedure;
for d in [2..5] do
 if d eq 5 then
  assert &and[Degree(c) lt 5:c in rawcs];AF:=AlgorithmicFunctionField(FunctionField(C));xs:=[CurvePlace(C,x):x in Places(AF,5)];assert #xs eq 74688;
  print "DECIC612_AF_TARGET_ENUM_COMPLETE",#xs,Cputime(tt);rr,labels,keys:=Decic612OrbitsCompact(X,xs);
 else xs:=NoncuspidalPlaces(X,d);assert #xs gt 0;rr,labels,keys:=Decic612Orbits(X,xs);end if;
 rcl:=[m2(Divisor(r)-d*Divisor(c0)):r in rr];Append(~allxs,xs);Append(~reps,rr);Append(~repc,rcl);
 if d le 3 then Append(~allcl,[m2(Divisor(x)-d*Divisor(c0)):x in xs]);end if;
 comp:=WeightedVectors(rawds,5-d);cc:=[&+[v[i]*rc[i]:i in [1..#rawcs]]:v in comp];
 for i in [1..#rr] do for j in [1..#comp] do Check(rcl[i]+cc[j],[*d,i,comp[j]*],~targets,~matches);end for;end for;
 Append(~counts,[d,#xs,#rr,#comp]);Append(~coverage,[*d,keys,labels,[Index(xs,r):r in rr]*]);
 print "DECIC612_D5_TARGET_PROGRESS",d,#xs,#rr,#targets,#matches,Cputime(tt);
 Write(Sprintf("work/deg10/continuation_20261029/diagonal612_D5_v4_target%o_checkpoint.m",d),Sprint([*targets,matches,counts,coverage*]));
end for;
comp:=WeightedVectors(rawds,1);cc:=[&+[v[i]*rc[i]:i in [1..#rawcs]]:v in comp];
for i in [1..#repc[1]] do
 for j in [1..#allcl[1]] do for k in [1..#comp] do Check(repc[1][i]+allcl[1][j]+cc[k],[*22,i,j,comp[k]*],~targets,~matches);end for;end for;
 for j in [1..#allcl[2]] do Check(repc[1][i]+allcl[2][j],[*23,i,j*],~targets,~matches);end for;
end for;
Write("work/deg10/continuation_20261029/diagonal_6_12_D5_arithmetic.txt",Sprint(prefix cat [*rawds,[Eltseq(c):c in rc],targets,matches,counts,coverage,[[Eltseq(c):c in cs]:cs in allcl]*]));
print "DECIC612_DIAGONAL_D5_COMPLETE",#source,#targets,#matches,counts,Cputime(tt);quit;
