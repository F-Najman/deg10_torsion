SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(7*10^9);Alarm(2100);tt:=Cputime();AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");Attach("work/deg10/continuation_20261031b/C48/inputs/fast_hecke.m");load "work/deg10/continuation_20261031b/C48/inputs/pair_af_place_linear.m";root:="work/deg10/continuation_20261031b/C48/";X:=MDX1(48,GF(5));F0,mp:=AlgorithmicFunctionField(FunctionField(Curve(X)));F:=RationalExtensionRepresentation(F0);assert Genus(F)eq 37;xy:=[F!mp(v):v in X`_coordinates_xy];bcf:=[F!mp(v):v in X`_coordinates];OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);k:=GF(5);R:=PolynomialRing(k);src:=eval Read(root cat "low_hecke_orbits.m");fs:=AssociativeArray();for r in src[1]do fs[r[1]]:=ext<k|R!r[2]>;end for;
function Get(E,P)
 vals:=TateNormalForm_xy(E,P);try pt:=C48PairAFPlaceLinear(xy,vals);catch err pt:=C48PairAFPlaceLinear(bcf,TateNormalForm_bc(E,P));end try;return pt;
end function;
function ER(v)return [[Integers()!a:a in Coefficients(Numerator(v))],[Integers()!a:a in Coefficients(Denominator(v))]];end function;
function EF(v)return [ER(z):z in Eltseq(F!v)];end function;
function PD(p)return [*Degree(p),Order(Ideal(p))eq OF select 1 else 0,[EF(z):z in Basis(Ideal(p))]*];end function;
Write(root cat TAG cat "_model.m",Sprint([*[ER(v):v in Coefficients(DefiningPolynomial(F))],[EF(v):v in xy],[EF(v):v in bcf]*]):Overwrite:=true);inputs:=[r:r in src[3]|(r[4]eq 0)eq PURE];out:=[* *];pointrows:=[* *];cache:=AssociativeArray();left:=[];print "D31B_C48_LOW_SOURCE37_READY",Q,SHARD,Degree(F),#inputs,Cputime(tt);
for ci in[SHARD+1..#inputs by NSHARDS]do cas:=inputs[ci];id:=cas[1];A:=DivisorGroup(F)!0;
 for ri in cas[3]do
  if not IsDefined(cache,ri)then
   rec:=src[2][ri];assert rec[1]eq ri;d:=rec[2];ky:=rec[3];K:=fs[d];b:=Evaluate(R!ky[1..d],K.1);c:=Evaluate(R!ky[d+1..2*d],K.1);E:=EllipticCurve([K|1-c,-b,-b,0,0]);P:=E![0,0];assert Order(P)eq 48;ps:=Get(E,P);pd:=Get(E,Q*P);assert d mod Degree(ps)eq 0 and d mod Degree(pd)eq 0;ns:=d div Degree(ps);nd:=Q*d div Degree(pd);B:=-ns*Divisor(ps)-nd*Divisor(pd);imgs:=[* *];kernelcount:=0;
   for z in MDIsogeniesFast(E,P,Q)do EE,PP,e:=Explode(z);p:=Get(EE,PP);assert d*e mod Degree(p)eq 0;n:=d*e div Degree(p);B+:=n*Divisor(p);Append(~imgs,[*e,n,PD(p)*]);kernelcount+:=e;end for;
   assert kernelcount eq Q+1 and Degree(B)eq 0;cache[ri]:=B;Append(~pointrows,[*ri,d,Q,PD(ps),ns,PD(pd),nd,imgs*]);print "D31B_C48_LOW_HEK_POINT",Q,ri,#pointrows,Cputime(tt);
  end if;
  A+:=cache[ri];
 end for;
 assert Degree(A)eq 0;hh:=#Basis(A);assert hh in[0,1];Append(~out,[*id,cas[2],cas[3],cas[4],hh*]);if hh eq 1 then Append(~left,id);end if;Write(root cat TAG cat "_progress.m",Sprint([*pointrows,out,left*]):Overwrite:=true);print "D31B_C48_LOW_HEK_CASE",Q,id,hh,#out,#left,Cputime(tt);
end for;
Write(root cat TAG cat ".m",Sprint([*pointrows,out,left*]):Overwrite:=true);print "D31B_C48_LOW_HEK_BATCH_COMPLETE",Q,SHARD,#out,#pointrows,left,Cputime(tt);quit;
