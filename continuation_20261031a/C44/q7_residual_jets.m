SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(12*10^9);tt:=Cputime();
AttachSpec("work/degree9/shared/mdmagma/v2/mdmagma.spec");Attach("work/degree9/shared/fast_hecke.m");
load "work/deg10/continuation_20261031a/C44/inputs/canonical_data.m";
assert AFBaseCoordinate eq 2;PF:=PolynomialRing(GF(3));X:=MDX1(44,GF(3));C:=Curve(X);cf:=DefiningPolynomial(C);AA<xx,yy>:=AffineSpace(GF(3),2);AR:=CoordinateRing(AA);
if Rank(Parent(cf)) eq 3 then fxy:=Evaluate(cf,[AR.1,AR.2,AR!1]);else fxy:=Evaluate(cf,[AR.1,AR.2]);end if;fx:=Derivative(fxy,1);fy:=Derivative(fxy,2);
function Normalize(z)
 k:=Index([w ne 0:w in z],true);assert k gt 0;return [w/z[k]:w in z];
end function;
function CanonicalValues(xy)
 LF:=Parent(xy[1]);tv:=xy[2];gv:=xy[1];gp:=[LF!1];for i in [1..18] do Append(~gp,gp[#gp]*gv);end for;zz:=[];
 for co in CanonicalCoordinateData do val:=LF!0;for j in [1..#co] do val+:=(Evaluate(PF!co[j][1],tv)/Evaluate(PF!co[j][2],tv))*gp[j];end for;Append(~zz,val);end for;return Normalize(zz);
end function;
function IntCode(v)
 z:=0;pw:=1;for c in v do z+:=(Integers()!c)*pw;pw*:=3;end for;return z;
end function;
function ReduceXY(xy)
 LF:=Parent(xy[1]);DD:=Degree(LF);dd:=Min([r:r in Divisors(DD)|&and[w^(3^r) eq w:w in xy]]);
 candidates:=[xy[1],xy[2],xy[1]+xy[2],xy[1]+2*xy[2]];ix:=0;
 for i in [1..4] do mp:=MinimalPolynomial(candidates[i],GF(3));if Degree(mp) eq dd then ix:=i;theta:=candidates[i];break;end if;end for;assert ix gt 0;
 BB:=Matrix(GF(3),[Eltseq(theta^i,GF(3)):i in [0..dd-1]]);assert Nrows(BB) eq dd and Ncols(BB) eq DD and Rank(BB) eq dd;
 cx:=Eltseq(Solution(BB,Vector(GF(3),Eltseq(xy[1],GF(3)))));cy:=Eltseq(Solution(BB,Vector(GF(3),Eltseq(xy[2],GF(3)))));
 assert Evaluate(PF!cx,theta) eq xy[1] and Evaluate(PF!cy,theta) eq xy[2];
 RF<ar>:=ext<GF(3)|mp>;xy:=[Evaluate(PF!cx,ar),Evaluate(PF!cy,ar)];assert Degree(RF) eq dd and Evaluate(fxy,xy) eq 0;
 z:=CanonicalValues(xy);key:=<dd,ix,IntCode(Eltseq(mp)),IntCode(cx),IntCode(cy)>;return xy,z,dd,key;
end function;
function EncodingBasis(LF)
 dd:=Degree(LF);theta:=LF.1;if Degree(MinimalPolynomial(theta,GF(3))) ne dd then theta:=PrimitiveElement(LF);end if;
 BB:=Matrix(GF(3),[Eltseq(theta^i,GF(3)):i in [0..dd-1]]);assert Nrows(BB) eq dd and Ncols(BB) eq dd and Rank(BB) eq dd;return [Integers()!c:c in Eltseq(MinimalPolynomial(theta,GF(3)))],BB^-1;
end function;
function Enc(w,IB)return [Integers()!a:a in Eltseq(Vector(GF(3),Eltseq(w,GF(3)))*IB)];end function;
function JetData(xy,zz,mm)
 LF:=Parent(xy[1]);prec:=Max(mm,2);PS<T>:=PowerSeriesRing(LF,prec);xs:=PS!xy[1];ys:=PS!xy[2];vx:=Evaluate(fx,xy);vy:=Evaluate(fy,xy);
 if vy ne 0 then dir:=1;xs+:=T;for k in [1..prec-1] do er:=Coefficient(Evaluate(fxy,[xs,ys]),k);ys-:=(er/vy)*T^k;end for;
 else assert vx ne 0;dir:=2;ys+:=T;for k in [1..prec-1] do er:=Coefficient(Evaluate(fxy,[xs,ys]),k);xs-:=(er/vx)*T^k;end for;end if;
 assert &and[Coefficient(Evaluate(fxy,[xs,ys]),k) eq 0:k in [0..prec-1]];
 gp:=[PS!1];for i in [1..18] do Append(~gp,gp[#gp]*xs);end for;zs:=[];
 for co in CanonicalCoordinateData do val:=PS!0;for j in [1..#co] do val+:=(Evaluate(PF!co[j][1],ys)/Evaluate(PF!co[j][2],ys))*gp[j];end for;Append(~zs,val);end for;
 ki:=Index([Coefficient(w,0) ne 0:w in zs],true);assert ki gt 0;zs:=[w/zs[ki]:w in zs];assert [Coefficient(w,0):w in zs] eq zz and exists{w:w in zs|Coefficient(w,1) ne 0};
 fp,IB:=EncodingBasis(LF);xycode:=[Enc(w,IB):w in xy];xyjets:=[[Enc(Coefficient(w,k),IB):k in [0..mm-1]]:w in [xs,ys]];zjets:=[[Enc(Coefficient(w,k),IB):k in [0..mm-1]]:w in zs];return fp,xycode,xyjets,zjets,dir;
end function;

SelectedData:=[* *];
Append(~SelectedData,[*[1, 2, 0, 2, 0, 0, 0, 1, 1, 2, 1],[[1, 1, 2, 1, 0, 2, 0, 2, 2, 0], [1, 1, 1, 0, 1, 2, 1, 2, 1, 1], [1, 1, 1, 0, 1, 2, 1, 2, 1, 1], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],[[0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]*]);
Append(~SelectedData,[*[2, 0, 0, 1, 1, 1, 0, 0, 2, 2, 1],[[1, 0, 1, 1, 2, 2, 2, 0, 0, 1], [0, 0, 0, 0, 1, 2, 2, 2, 0, 0], [0, 0, 0, 0, 1, 2, 2, 2, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],[[0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]*]);

cases:=[* *];rf:=recformat<P:PtEll>;
for row in SelectedData do
 L<ar>:=ext<GF(3)|PF!row[1]>;EE:=EllipticCurve([Evaluate(PF!a,ar):a in row[2]]);PP:=EE![Evaluate(PF!a,ar):a in row[3]];
 assert 44*PP eq EE!0 and 22*PP ne EE!0 and 4*PP ne EE!0;
 lp:=rec<rf|P:=PP>;Q:=ModuliPoint(X,EE,lp);assert Degree(Q) eq 10;Append(~cases,[Q]);
end for;
units:=Sort(Setseq({Min(a,44-a):a in [1..43]|GCD(a,44) eq 1}));
out:=Open("work/deg10/continuation_20261031a/C44/q7_finite_jets.raw.jsonseq","w");
fprintf out,"{\"type\":\"header\",\"N\":44,\"p\":3,\"q\":7,\"original_cases\":[510,511],\"requested_cases\":2}\n";
sourceKeys10:={};sourceKeys5:={};
for k in [1..#cases] do
 PPcases:=cases[k];total:=&+[Degree(P):P in PPcases];tg:=AssociativeArray();ug:=AssociativeArray();isoSizes:=[* *];
 for si in [1..#PPcases] do
  Pplace:=PPcases[si];ds:=Degree(Pplace);E:=EllipticCurve(X,Pplace);PP:=LevelStructure(X,Pplace)`P;assert 44*PP eq E!0 and 22*PP ne E!0 and 4*PP ne E!0;
  K0:=BaseRing(E);assert Degree(K0) eq ds;fp0,IB0:=EncodingBasis(K0);
  fprintf out,"{\"type\":\"source\",\"case\":%o,\"source_index\":%o,\"source_degree\":%o,\"field_polynomial\":%o,\"E_ai\":%o,\"P_xy\":%o}\n",k,si,ds,fp0,[Enc(a,IB0):a in aInvariants(E)],[Enc(PP[i],IB0):i in [1,2]];
  for a in units do
   xy0,z0,d0,key0:=ReduceXY(TateNormalForm_xy(E,a*PP));assert d0 eq ds;
   if ds eq 10 then Include(~sourceKeys10,key0);else Include(~sourceKeys5,key0);end if;
   if a in [1,7] then
    mm:=a eq 1 select 1 else 7;rr:=[*xy0,z0,d0,mm,key0*];
    if IsDefined(ug,key0) then old:=ug[key0];old[4]+:=mm;ug[key0]:=old;else ug[key0]:=rr;end if;
   end if;
  end for;
  es:=[];for iso in MDIsogeniesFast(E,PP,7) do
   EG,PG,ee:=Explode(iso);assert 44*PG eq EG!0 and 22*PG ne EG!0 and 4*PG ne EG!0;
   xy0,z0,d0,key0:=ReduceXY(TateNormalForm_xy(EG,PG));assert (ee*ds) mod d0 eq 0;mm:=ee*ds div d0;rr:=[*xy0,z0,d0,mm,key0*];Append(~es,ee);
   if IsDefined(tg,key0) then old:=tg[key0];old[4]+:=mm;tg[key0]:=old;else tg[key0]:=rr;end if;
  end for;assert &+es eq 8;Append(~isoSizes,es);
 end for;
 assert &+[v[3]*v[4]:v in Values(tg)] eq 8*total and &+[v[3]*v[4]:v in Values(ug)] eq 8*total;
 fprintf out,"{\"type\":\"case\",\"case\":%o,\"support_degrees\":%o,\"isogeny_orbit_sizes\":%o,\"positive_degree\":%o,\"negative_degree\":%o}\n",k,[Degree(P):P in PPcases],isoSizes,8*total,8*total;
 for side in [1,2] do
  dict:=side eq 1 select tg else ug;label:=side eq 1 select "T" else "U";kk:=0;
  for key in Sort(Setseq(Keys(dict))) do
   kk+:=1;rr:=dict[key];fp,xycode,xyjets,zjets,dir:=JetData(rr[1],rr[2],rr[4]);
   fprintf out,"{\"type\":\"jet\",\"case\":%o,\"side\":\"%o\",\"support\":%o,\"degree\":%o,\"multiplicity\":%o,\"field_polynomial\":%o,\"xy\":%o,\"xy_jets\":%o,\"canonical_jets\":%o,\"uniformizer_coordinate\":%o}\n",k,label,kk,rr[3],rr[4],fp,xycode,xyjets,zjets,dir;
  end for;
 end for;Flush(out);
 if k mod 25 eq 0 then printf "D31A_C44_Q7_JET_CASE index=%o/%o T_supports=%o U_supports=%o time=%o\n",k,#cases,#Keys(tg),#Keys(ug),Cputime(tt);end if;
end for;
assert #sourceKeys5 eq 0;
fprintf out,"{\"type\":\"complete\",\"cases\":%o,\"distinct_source10\":%o,\"distinct_source5\":%o}\n",#cases,#sourceKeys10,#sourceKeys5;delete out;
printf "D31A_C44_Q7_FINITE_JETS_COMPLETE cases=%o distinct10=%o distinct5=%o time=%o\n",#cases,#sourceKeys10,#sourceKeys5,Cputime(tt);quit;
