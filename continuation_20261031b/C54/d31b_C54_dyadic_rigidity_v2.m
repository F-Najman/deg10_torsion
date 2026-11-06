SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(6*10^9);Alarm(1100);tt:=Cputime();AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");Attach("work/deg10/continuation_20261031b/C54/inputs/class_rational_reduction.m");Attach("work/deg10/continuation_20261031b/C54/inputs/fast_hecke.m");root:="work/deg10/continuation_20261031b/C54/";X:=MDX1(27,GF(2));C:=Curve(X);FF:=FunctionField(C);AF,am:=AlgorithmicFunctionField(FF);F:=RationalExtensionRepresentation(AF);assert Genus(F)eq 13;cs:=Cusps(X);low:=[p:p in cs|Degree(p)le 4];sig:=[CuspSignature(X,p):p in low];assert #Seqset(sig)eq #low;
EE:=EllipticCurve([FF!z:z in X`_E]);P:=EE![FF!z:z in X`_P];bc:=TateNormalForm_bc(EE,2*P);b,c:=Explode(bc);units:=[b/(16*b^2+(1-20*c-8*c^2)*b+c*(c-1)^3),b];perm:=[];for p in low do sg:=[Valuation(v,p):v in units];ix:=Index(sig,sg);assert ix gt 0 and Degree(low[ix])eq Degree(p);Append(~perm,ix);end for;assert #Seqset(perm)eq #low;
c0:=[p:p in low|Degree(p)eq 1][1];A:=FunctionFieldDivisor(Divisor(c0));G,_,phi:=ClassGroup(F:ReductionDivisor:=A,Proof:=true);assert Invariants(G)eq[3,3,52497,0];cp:=[Support(FunctionFieldDivisor(Divisor(p)))[1]:p in low];cd:=[Degree(p):p in cp];cg:=[phi(Divisor(p)-Degree(p)*A):p in cp];vectors:=AssociativeArray();for d in[0..4]do vectors[d]:=[];end for;
procedure Vectors(j,left,v,~out)
 if j gt #cd then if left eq 0 then Append(~out,v);end if;return;end if;
 for n in[0..left div cd[j]]do Vectors(j+1,left-n*cd[j],v cat[n],~out);end for;
end procedure;
for d in[0..4]do vv:=[];Vectors(1,d,[],~vv);vectors[d]:=vv;end for;
rows:=[* *];left:=[* *];ncinfo:=[* *];tot:=0;kill:=0;
for d in[6..10]do ps:=NoncuspidalPlaces(X,d);ps:=[p:p in ps|jInvariant(EllipticCurve(X,p))ne 0];reps:=PlacesUpToDiamondFast(X,ps);print "D31B_C54_DYADIC_ORBITS",d,#ps,#reps,Cputime(tt);
 for ni in[1..#reps]do p:=Support(FunctionFieldDivisor(Divisor(reps[ni])))[1];N:=Divisor(p);ng:=phi(N-d*A);cache:=AssociativeArray();localrows:=[* *];e:=10-d;
  for ceDeg in[0..e]do for cv in vectors[ceDeg],av in vectors[e-ceDeg]do rv:=[av[j]+cv[j]:j in[1..#cd]];for j in[1..#cd]do rv[perm[j]]+:=cv[j];end for;
   if IsDefined(cache,rv)then h:=cache[rv];else D:=N+&+[rv[j]*Divisor(cp[j]):j in[1..#cd]];assert Degree(D)eq 10+ceDeg;h:=Dimension(D);cache[rv]:=h;end if;
   cl:=ng+&+[(av[j]+cv[j])*cg[j]:j in[1..#cd]];cl2:=ng+&+[av[j]*cg[j]+cv[j]*cg[perm[j]]:j in[1..#cd]];tot+:=1;row:=[*tot,d,ni,av,cv,h,Eltseq(cl),Eltseq(cl2)*];Append(~localrows,row);if h eq 1 then kill+:=1;else Append(~left,row);end if;
  end for;end for;rp:=RepresentativePoint(reps[ni]);kk:=Parent(rp[1]);assert #kk eq 2^d;Append(~ncinfo,[*d,ni,[Integers()!z:z in Coefficients(DefiningPolynomial(kk,GF(2)))],[[Integers()!z:z in Eltseq(v,GF(2))]:v in Eltseq(rp)],#cache*]);rows cat:=localrows;print "D31B_C54_DYADIC_NC_CASE",d,ni,#localrows,#cache,#[v:v in localrows|v[6]ne 1],Cputime(tt);Write(root cat "dyadic_rigidity_v2_progress.m",Sprint([*Invariants(G),[Eltseq(v):v in Rows(RelationMatrix(G))],sig,cd,perm,[Eltseq(v):v in cg],rows,left,ncinfo*]):Overwrite:=true);
 end for;
end for;
Write(root cat "dyadic_rigidity_v2.m",Sprint([*Invariants(G),[Eltseq(v):v in Rows(RelationMatrix(G))],sig,cd,perm,[Eltseq(v):v in cg],rows,left,ncinfo*]):Overwrite:=true);print "D31B_C54_DYADIC_RIGIDITY_COMPLETE",tot,kill,#left,Cputime(tt);quit;
