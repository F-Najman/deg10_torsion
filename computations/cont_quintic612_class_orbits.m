SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(8*10^9);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
Attach("work/deg10/continuation_20261029/deps/fast_hecke.m");
Attach("work/deg10/continuation_20261029/deps/fast_hecke_11.m");
Attach("work/deg10/continuation_20261029/deps/fast_diamond_hash_m2.m");
Attach("work/deg10/continuation_20261029/deps/model_6_12.m");
Attach("work/deg10/continuation_20261029/deps/full_aut_6_12.m");
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
data:=eval Read("work/deg10/continuation_20261029/diagonal612_D5_source_checkpoint_v5.m");ll:=eval Read("work/deg10/continuation_20261029/M6N12_D5_all_source_lifts.m");S:=eval Read("work/deg10/continuation_20261029/RR612_norm19_survivor_indices.txt");
FB:=FreeAbelianGroup(#data[5][1][1]);G,pi:=quo<FB|[FB!r:r in data[4][1]]>;cg:=[pi(FB!v):v in data[5][1]];
keymap:=AssociativeArray();for k in S do v:=ll[3][k][2];cl:=&+[v[j]*cg[j]:j in [1..24]];assert not IsDefined(keymap,cl);keymap[cl]:=k;end for;
point:=AssociativeArray();for sh in [0..7] do
 cert:=eval Read(Sprintf("work/deg10/continuation_20261029/cont_source612_recert_v2_shard%o_arithmetic.txt",sh));for row in cert do if row[1] eq 5 then point[row[2]]:=pi(FB!row[3]);end if;end for;
end for;assert #Keys(point) eq 3264;
X:=MDX612(GF(7));xs:=NoncuspidalPlaces(X,5);assert #xs eq 3264;loc:=AssociativeArray();for i in [1..#xs] do loc[xs[i]]:=i;end for;
remaining:=Seqset([k:k in S|k lt 769 or k gt 864]);assert #remaining eq 96;orbits:=[* *];pars:=Decic612Parameters();
while #remaining gt 0 do
 k:=Minimum(remaining);ids:=ll[3][k][1];sid:=[i:i in ids|data[8][i][1] eq 5][1];i:=data[8][sid][2];assert 2*point[i] eq &+[ll[3][k][2][j]*cg[j]:j in [1..24]];
 x:=xs[i];E:=EllipticCurve(X,x);ls:=LevelStructure612(X,x);P:=ls`P;Q:=ls`Q;images:=[];classes:=[];keys:=[];K:=BaseRing(E);
 for par in pars do
  a,b,c,d:=Explode(par);PP:=a*P+2*b*Q;QQ:=c*P+d*Q;assert Order(PP) eq 6 and Order(QQ) eq 12 and WeilPairing(PP,2*QQ,6) eq WeilPairing(P,2*Q,6);
  level:=rec<recformat<P:PtEll,Q:PtEll>|P:=PP,Q:=QQ>;im:=ModuliPoint612(X,E,level);assert IsDefined(loc,im);j:=loc[im];cl:=2*point[j];assert IsDefined(keymap,cl);kk:=keymap[cl];
  Append(~images,j);Append(~classes,kk);Append(~keys,MDM2FrobeniusKey(Decic612Normalize(E,PP,QQ),K,7,5));
 end for;
 orb:=Seqset(classes);assert orb subset remaining;remaining diff:=orb;Append(~orbits,[*k,i,images,classes,keys*]);print "DECIC612_QUINTIC_CLASS_ORBIT",k,#orb,Cputime(tt);
end while;
Write("work/deg10/continuation_20261029/M6N12_quintic_class_orbits.m",Sprint([*pars,orbits*]));print "DECIC612_QUINTIC_CLASS_ORBITS_COMPLETE",#orbits,Cputime(tt);quit;
