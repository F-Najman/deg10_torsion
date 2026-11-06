/* Fresh octic M4 two-multiple diagonal. Global C is the subgroup of
   F-rational cusp DIVISORS, not assumed to be the full Jacobian torsion. */
SetColumns(0);SetNthreads(1);SetQuitOnError(true);SetMemoryLimit(8*10^9);t0:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
Attach("work/deg10/continuation_20261029/deps/fast_hecke.m");Attach("work/deg10/continuation_20261029/deps/fast_hecke_11.m");
Attach("work/deg10/continuation_20261029/deps/fast_diamond_hash_m2.m");
Attach("work/deg10/continuation_20261029/deps/fast_full_aut_hash_m4_16.m");
Attach("work/deg10/continuation_20261029/deps/fast_full_aut_hash_m4_16_relative.m");
Attach("work/deg10/shared/deps/class_rational_reduction.m");
function CuspUnitsM4(X)
 FF:=FunctionField(Curve(X));E:=EllipticCurve([FF!f:f in X`_E]);P:=E![FF!f:f in X`_P];Q:=E![FF!f:f in X`_Q];
 ai:=aInvariants(E);B:=ai[2];assert ai eq [FF!1,B,B,FF!0,FF!0] and P eq E![0,0];
 assert Discriminant(E) eq B^4*(1-16*B);fs:=[B,1-16*B];
 for k in [1..8] do for j in [0..3] do
  A:=k*Q+j*P;assert A ne E!0 and A[1] ne 0 and A[2] ne 0 and A[1]+B ne 0;
  fs cat:=[A[1],A[2],A[1]+B,A[2]+A[1]+B];
 end for;end for;
 for k in [1..8] do A:=k*Q;for j in [1..3] do R:=A+j*P;assert A ne R and A ne -R;Append(~fs,A[1]-R[1]);end for;end for;
 rawcs:=Cusps(X);groups:=AssociativeArray();
 for c in rawcs do sig:=[Valuation(f,c):f in fs];
  if IsDefined(groups,sig) then Append(~groups[sig],c);else groups[sig]:=[c];end if;
 end for;
 sigs:=Sort(Setseq(Keys(groups)));assert #sigs eq 24;
 return [&+[Divisor(c):c in groups[sig]]:sig in sigs],sigs;
end function;
Free:=FreeAbelianGroup(24);Xs:=[* *];Gs:=[* *];hs:=[* *];ms:=[* *];css:=[* *];bases:=[* *];cgs:=[* *];fps:=[];
sigs:=[];ds:=[];
for p in [5,9] do
 X:=MDX11(4,16,GF(p));C:=Curve(X);assert Genus(X) eq 13;cs,ss:=CuspUnitsM4(X);dd:=[Degree(c):c in cs];
 assert {*d:d in dd*} eq {*1^^16,2^^4,4^^4*};
 if p eq 5 then sigs:=ss;ds:=dd;else assert ss eq sigs and dd eq ds;end if;
 c0:=Support(cs[Index(ds,1)])[1];print "M4N16_D4_DIAGONAL_CUSP_SETUP",p,ds,Cputime(t0);
 G,m1,m2:=ClassGroupRationalReduction(C,c0:Proof:=true);
 J:=JH(64,[17]);assert Dimension(J) eq 13;fp:=FrobeniusPolynomial(J,p eq 9 select 3 else 5);expected:=p eq 9 select Evaluate(fp,1)*Evaluate(fp,-1) else Evaluate(fp,1);assert &*[n:n in Invariants(G)|n ne 0] eq Integers()!expected;
 hg:=[m2(cs[i]-ds[i]*Divisor(c0)):i in [1..24]];h:=hom<Free->G|hg>;H:=Image(h);
 assert Invariants(H) eq [2,2,2,4,4,20,40,40,80];
 if p eq 5 then assert Invariants(quo<G|H>) eq [2,2,2,2,2,0];end if;
 Append(~Xs,X);Append(~Gs,G);Append(~hs,h);Append(~ms,m2);Append(~css,cs);Append(~bases,c0);Append(~cgs,hg);Append(~fps,Coefficients(fp));
 print "M4N16_D4_DIAGONAL_GROUP",p,Invariants(G),Invariants(H),Invariants(quo<G|H>),Cputime(t0);
end for;
h5:=hs[1];h9:=hs[2];G5:=Gs[1];G9:=Gs[2];
for v in Generators(Kernel(h5)) do assert h9(Free!v) eq G9!0;end for;
for v in Generators(Kernel(h9)) do assert h5(Free!v) eq G5!0;end for;

function WeightedVectors(ds,d)
 if #ds eq 0 then return d eq 0 select [[]] else [];end if;
 out:=[];for k in [0..d div ds[1]] do for v in WeightedVectors(ds[2..#ds],d-k*ds[1]) do Append(~out,[k] cat v);end for;end for;return out;
end function;
X:=Xs[1];C:=Curve(X);c0:=bases[1];m2:=ms[1];source:=[* *];cache:=AssociativeArray();predictions:=AssociativeArray();sourcecoverage:=[* *];
for d in [3..5] do
 xs:=NoncuspidalPlaces(X,d);assert #xs eq [64,160,512][d-2];
 complement:=WeightedVectors(ds,5-d);cc:=[&+[v[i]*cgs[1][i]:i in [1..24]]:v in complement];
 for i in [1..#xs] do
  a:=m2(Divisor(xs[i])-d*Divisor(c0));
  for j in [1..#complement] do
   cl:=a+cc[j];tw:=2*cl;assert tw in Image(h5);key:=Eltseq(tw);
   if not IsDefined(cache,key) then cache[key]:=h9(Free!(tw@@h5));end if;
   pred:=cache[key];Append(~source,[*d,i,complement[j],Eltseq(cl),Eltseq(pred)*]);ky:=Eltseq(pred);
   if IsDefined(predictions,ky) then Append(~predictions[ky],#source);else predictions[ky]:=[#source];end if;
  end for;
 end for;
 Append(~sourcecoverage,[*d,#xs,#complement*]);print "DECIC416_DIAGONAL_SOURCE_PROGRESS",d,#source,#Keys(predictions),Cputime(t0);
end for;
assert #source eq 12032;print "DECIC416_DIAGONAL_SOURCE_COMPLETE",#source,#Keys(predictions),Cputime(t0);
X:=Xs[2];C:=Curve(X);c0:=bases[2];m2:=ms[2];rawcs:=Cusps(X);rawds:=[Degree(c):c in rawcs];
assert {*d:d in rawds*} eq {*1^^24,2^^8*};
rc:=[m2(Divisor(c)-Degree(c)*Divisor(c0)):c in rawcs];
targets:=[* *];matches:=[* *];coverage:=[* *];counts:=[];
for d in [3..5] do
 xs:=NoncuspidalPlaces(X,d);rr,labels,keys:=PlacesUpToFullAutM4SixteenRelative(X,xs);
 complement:=WeightedVectors(rawds,5-d);cc:=[&+[v[i]*rc[i]:i in [1..#rawcs]]:v in complement];
 for i in [1..#rr] do
  a:=m2(Divisor(rr[i])-d*Divisor(c0));
  for j in [1..#complement] do
   cl:=a+cc[j];Append(~targets,[*d,i,complement[j],Eltseq(cl)*]);ky:=Eltseq(2*cl);
   if IsDefined(predictions,ky) then Append(~matches,[*#targets,predictions[ky]*]);end if;
  end for;
 end for;
 Append(~coverage,[*d,keys,labels,[Index(xs,x):x in rr]*]);Append(~counts,[d,#xs,#rr,#complement]);
 print "DECIC416_DIAGONAL_TARGET_PROGRESS",d,#xs,#rr,#targets,#matches,Cputime(t0);
end for;
Write("work/deg10/continuation_20261029/diagonal_4_16_D5_arithmetic.txt",Sprint([*sigs,ds,Invariants(G5),[Eltseq(r):r in Rows(RelationMatrix(G5))],Invariants(G9),[Eltseq(r):r in Rows(RelationMatrix(G9))],[[Eltseq(v):v in cg]:cg in cgs],source,targets,matches,counts,coverage,[Eltseq(a):a in rc],fps,sourcecoverage*]));
print "DECIC416_DIAGONAL_D5_COMPLETE",#source,#targets,#matches,counts,Cputime(t0);quit;
