/* Necessary three-multiple diagonal test. ALL F4 effective degree4 divisors
   are included. AtF7 the entirely cuspidal branch has a global F-cusp lift
   and is excluded separately using odd-prime injectivity and gon_F>=6.
   Noncuspidal F7 patterns are single2/3/4 with every cusp complement,
   and2+2 including repetition. No global full-generation claim is used. */
SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(4*10^9);t0:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
Attach("work/deg10/shared/deps/class_rational_reduction.m");
Free:=FreeAbelianGroup(20);Xs:=[* *];Gs:=[* *];Hs:=[* *];hs:=[* *];ms:=[* *];
css:=[* *];basecs:=[* *];signatures:=[];degrees:=[];
function WeightedVectors(ds,d)
 if #ds eq 0 then return d eq 0 select [[]] else [];end if;
 result:=[];
 for k in [0..d div ds[1]] do
  for v in WeightedVectors(ds[2..#ds],d-k*ds[1]) do Append(~result,[k] cat v);end for;
 end for;
 return result;
end function;
for p in [4] do
 X:=MDX11(3,15,GF(p));C:=Curve(X);FF:=FunctionField(C);
 E:=EllipticCurve([FF!f:f in X`_E]);P:=E![FF!f:f in X`_P];Q:=E![FF!f:f in X`_Q];
 ai:=aInvariants(E);assert ai[2] eq 0 and ai[4] eq 0 and ai[5] eq 0;
 assert P eq E![0,0];fs:=[ai[3],ai[1]^3-27*ai[3]];
 for k in [1..7] do
  for j in [0..2] do
   A:=k*Q+j*P;assert A[3] ne 0 and A[1] ne 0 and A[2] ne 0;
   Append(~fs,A[1]);Append(~fs,A[2]);
  end for;
  A:=k*Q;B:=A+P;assert A ne B and A ne -B;Append(~fs,A[1]-B[1]);
 end for;
 groups:=AssociativeArray();cs:=Cusps(X);
 for cusp in cs do
  sig:=[Valuation(f,cusp):f in fs];
  if IsDefined(groups,sig) then Append(~groups[sig],cusp);else groups[sig]:=[cusp];end if;
 end for;
 sigs:=Sort(Setseq(Keys(groups)));assert #sigs eq 20;
 blocks:=[&+[Divisor(c):c in groups[sig]]:sig in sigs];
 ds:=[Degree(o):o in blocks];assert {*d:d in ds*} eq {*1^^16,4^^4*};
 if p eq 4 then signatures:=sigs;degrees:=ds;
 else assert sigs eq signatures and ds eq degrees;end if;
 c0:=Support(blocks[Index(ds,1)])[1];
 print "DECIC315_DIAGONAL_CUSP_BLOCKS",p,#sigs,ds,Cputime(t0);
 G,m1,m2:=ClassGroupRationalReduction(C,c0:Proof:=true);
 gs:=[m2(blocks[i]-ds[i]*Divisor(c0)):i in [1..20]];
 h:=hom<Free->G|gs>;H:=Image(h);
 assert Invariants(H) eq [3,3,24,96,96];
 J:=JH(45,[16]);assert Dimension(J) eq 9;
 fp:=FrobeniusPolynomial(J,p eq 4 select 2 else 7);
 expected:=(p eq 4) select Evaluate(fp,1)*Evaluate(fp,-1) else Evaluate(fp,1);
 assert &*[n:n in Invariants(G)|n ne 0] eq Integers()!expected;
 Append(~Xs,X);Append(~Gs,G);Append(~Hs,H);Append(~hs,h);Append(~ms,m2);
 Append(~css,cs);Append(~basecs,c0);
 print "DECIC315_DIAGONAL_GROUP",p,Invariants(G),Invariants(H),Cputime(t0);
end for;


G:=Gs[1];m2:=ms[1];c0:=basecs[1];cs:=css[1];h:=hs[1];cl:=[m2(Divisor(c)-Degree(c)*Divisor(c0)):c in cs];
V:=VectorSpace(GF(3),6);vs:=[V!Eltseq(a)[1..6]:a in cl];perm:=[0:i in cs];
for key in Keys(groups) do
 inds:=[Index(cs,c):c in groups[key]];
 if #inds eq 1 then perm[inds[1]]:=inds[1];else assert #inds eq 2 and &and[Degree(cs[i]) eq 2:i in inds];perm[inds[1]]:=inds[2];perm[inds[2]]:=inds[1];end if;
end for;
assert Sort(perm) eq [1..24];inds:=[];W:=sub<V|>;
for i in [1..24] do if vs[i] notin W then Append(~inds,i);W:=sub<V|Basis(W) cat [vs[i]]>;end if;end for;
print "DECIC315_LOCAL_CUSP_THREE_RANK",Dimension(W);
if Dimension(W) eq 6 then
 A:=Matrix([vs[i]:i in inds]);B:=Matrix([vs[perm[i]]:i in inds]);M:=A^-1*B;assert &and[vs[i]*M eq vs[perm[i]]:i in [1..24]] and M^2 eq IdentityMatrix(GF(3),6);
 globalvs:=[V!Eltseq(h(Free.i))[1..6]:i in [1..20]];WG:=sub<V|globalvs>;WI:=Nullspace(M-IdentityMatrix(GF(3),6));
 print "DECIC315_GALOIS_FIXED_THREE",Dimension(WI),Dimension(WG),WI eq WG;
 Write("work/deg10/exclusions/315d5/galois3_data.txt",Sprint([*Invariants(G),[Eltseq(v):v in Rows(RelationMatrix(G))],[Degree(c):c in cs],perm,[Eltseq(v):v in vs],[Eltseq(v):v in globalvs],[Integers()!a:a in Eltseq(M)]*]));
end if;
print "DECIC315_GALOIS_THREE_DIAGNOSTIC_COMPLETE";quit;
