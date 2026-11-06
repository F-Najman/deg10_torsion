SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(3*10^9);tt:=Cputime();
ad:=eval Read("work/deg10/C236_resolution_20261030/enumeration/target18_aut_actions.m");Z:=FreeAbelianGroup(#ad[2][1]);Q,q:=quo<Z|[Z!r:r in ad[2]]>;assert Invariants(Q) eq [2,6];
function RawClass(v) assert #v eq 4 and v[4] eq 0;return &+[v[j]*(Q!ad[3][j]):j in [1..3]];end function;
autos:=[*hom<Q->Q|[Q!v:v in mm]>:mm in ad[4]*];for i in [1..Ngens(Q)] do assert autos[1](autos[2](Q.i)) eq autos[2](autos[1](Q.i));assert autos[3](autos[2](Q.i)) eq autos[2](autos[3](Q.i));end for;
function Power(f,e,v) for i in [1..e] do v:=f(v);end for;return v;end function;
function Key(vs,D) keys:=[];for i in [1..D] do Append(~keys,&cat[[Integers()!a:a in Eltseq(v,GF(5))]:v in vs]);vs:=[v^5:v in vs];end for;return Sort(keys)[1];end function;
function BC(E,P,Q) ai:=aInvariants(E);xx:=Q[1]-P[1];yy:=Q[2];return [xx^3/yy^2,(ai[2]+3*P[1])*xx^2/yy^2];end function;
for d in [3..7] do
 a:=eval Read(Sprintf("work/deg10/C236_resolution_20261030/enumeration/filter18_D%o_ownsmall.m",d));k:=GF(5);PP:=PolynomialRing(k);K<z>:=ext<k|PP!a[1]>;dict:=AssociativeArray();own:=[* *];
 for row in a[2] do
  ky:=row[3];b:=Evaluate(PP!ky[1..d],z);c:=Evaluate(PP!ky[d+1..2*d],z);E:=EllipticCurve([K|0,c,0,b*(1-b-c),0]);P:=E![0,0];S:=E![b,b];v1:=RawClass(row[4][1]);v2:=RawClass(row[4][2]);Append(~own,[*d,row[2],ky,[Eltseq(v1),Eltseq(v2)]*]);ims:={};
  for bb,cc in [0,1],u in [1,5,7,11,13,17] do
   e:=u in [1,17] select 0 else u in [5,13] select 1 else 2;w1:=Power(autos[1],bb,Power(autos[2],e,v1));w2:=Power(autos[3],cc,Power(autos[2],e,v2));im:=Key(BC(E,P+bb*(18*S),cc*P+u*S),d);Include(~ims,im);vals:=[Eltseq(w1),Eltseq(w2)];
   if IsDefined(dict,im) then assert dict[im] eq vals;else dict[im]:=vals;end if;
  end for;
  assert ky eq Sort(Setseq(ims))[1];
 end for;
 keys:=Sort(Setseq(Keys(dict)));assert #keys eq [48,120,720,2904,10416][d-2];small:=[*[*d,i,keys[i],dict[keys[i]]*]:i in [1..#keys]*];
 if d le 6 then
  old:=eval Read(Sprintf("work/deg10/C236_resolution_20261030/enumeration/filter18_D%o.m",d));assert #old[3] eq #small;
  for r in old[3] do assert [Eltseq(RawClass(v)):v in r[4]] eq dict[r[3]];end for;
 end if;
 Write(Sprintf("work/deg10/C236_resolution_20261030/enumeration/small_qclasses%o.m",d),Sprint([*a[1],own,small*]):Overwrite:=true);print "C236_SMALL_AUT_CLASS_PROPAGATION_COMPLETE",d,#own,#small,Cputime(tt);
end for;quit;
