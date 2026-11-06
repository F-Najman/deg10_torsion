// Root-allocated elliptic Hecke/jet stage; no function-field normalization.
SetColumns(0);SetNthreads(1);SetQuitOnError(true);SetMemoryLimit(2*10^9);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");Attach("work/deg10/continuation_20261030c/C65/inputs/fast_hecke.m");
load "work/deg10/continuation_20261030c/C65/inputs/canonical_unit_data.m";
function Code(z)return [Integers()!a:a in Eltseq(z,GF(2))];end function;
function SimplifyPoint(E,P)
 F:=BaseRing(E);a1,a2,a3,a4,a6:=Explode(aInvariants(E));assert a1 eq 1 and a3 eq 0;EE:=EllipticCurve([F!1,a2,F!0,F!0,a6+a4^2]);PP:=EE![P[1],P[2]+a4];assert jInvariant(EE)eq jInvariant(E)and 65*PP eq EE!0;assert TateNormalForm_xy(EE,PP)eq TateNormalForm_xy(E,P);return EE,PP;
end function;
function UnitValues(aa,xx,yy)
 a1,a2,a3,a4,a6:=Explode(aa);aa3:=2*yy+a3+a1*xx;aa2:=3*xx+a2;aa4:=3*xx^2+2*xx*a2+a4-a1*yy;aaa1:=2*aa4/aa3+a1;aaa2:=aa2-(aa4/aa3)^2-a1*aa4/aa3;b:=-aaa2^3/aa3^2;c:=-(aaa1*aaa2-aa3)/aa3;r:=b/c;s:=c^2/(b-c);t:=r*s-2*r+1;return[(s-r)/t,t/(s^2-s-r+1)];
end function;
function UnitMatrix(F,x0,y0)
 return VerticalJoin(Matrix(F,[[F!v[i]+x0*v[121+i]:i in[1..121]]:v in Rows(UnitXRelations)]),Matrix(F,[[F!v[i]+y0*v[121+i]:i in[1..121]]:v in Rows(UnitYRelations)]));
end function;
function CanonicalPoint(E,P)
 F:=BaseRing(E);xy:=TateNormalForm_xy(E,P);RR:=UnitMatrix(F,xy[1],xy[2]);assert Rank(RR)eq 120;zz:=Eltseq(Basis(Kernel(Transpose(RR)))[1]);ki:=Index([z ne 0:z in zz],true);return[z/zz[ki]:z in zz],xy;
end function;
function Psi65(xx,aa6)
 PS:=Parent(xx);pp:=[PS!0,PS!1,xx,xx^4+xx^3+aa6,xx^6+aa6*xx^2];
 for n in[5..65]do m:=n div 2;if IsOdd(n)then pn:=pp[m+3]*pp[m+1]^3+pp[m]*pp[m+2]^3;else pn:=pp[m+1]/xx*(pp[m+3]*pp[m]^2+pp[m-1]*pp[m+2]^2);end if;Append(~pp,pn);end for;return pp[66];
end function;
function CanonicalJets(E,P,z0)
 F:=BaseRing(E);aa:=aInvariants(E);assert aa[1]eq 1 and aa[3]eq 0 and aa[4]eq 0;PS<u>:=PowerSeriesRing(F,5);xx:=PS!P[1];yy:=PS!P[2];a2:=aa[2];a6:=PS!aa[5]+u;der:=Coefficient(Psi65(xx+u,PS!aa[5]),1);assert der ne 0 and P[1]ne 0;
 for n in[1..4]do xx+:=(Coefficient(Psi65(xx,a6),n)/der)*u^n;assert Valuation(Psi65(xx,a6))ge n+1;end for;
 for n in[1..4]do er:=Coefficient(yy^2+xx*yy+xx^3+a2*xx^2+a6,n);yy+:=(er/P[1])*u^n;assert Valuation(yy^2+xx*yy+xx^3+a2*xx^2+a6)ge n+1;end for;
 xy:=UnitValues([PS!1,PS!a2,PS!0,PS!0,a6],xx,yy);RM:=[UnitMatrix(F,Coefficient(xy[1],0),Coefficient(xy[2],0))];
 for n in[1..4]do Append(~RM,VerticalJoin(Matrix(F,[[Coefficient(xy[1],n)*v[121+i]:i in[1..121]]:v in Rows(UnitXRelations)]),Matrix(F,[[Coefficient(xy[2],n)*v[121+i]:i in[1..121]]:v in Rows(UnitYRelations)])));end for;
 ki:=Index([z ne 0:z in z0],true);assert z0[ki]eq 1;uk:=[i:i in[1..121]|i ne ki];M0:=Matrix(F,[[RM[1][r,c]:c in uk]:r in[1..Nrows(RM[1])]]);assert Rank(M0)eq 120;LL:=Solution(M0,IdentityMatrix(F,120));assert LL*M0 eq IdentityMatrix(F,120);zv:=[Vector(F,z0)];
 for n in[1..4]do
  br:=Vector(F,[0:i in[1..Nrows(M0)]]);for j in[1..n]do br+:=Vector(F,[&+[RM[j+1][r,c]*zv[n-j+1][c]:c in[1..121]]:r in[1..Nrows(M0)]]);end for;ur:=br*Transpose(LL);zz:=Vector(F,[0:i in[1..121]]);for j in[1..120]do zz[uk[j]]:=ur[j];end for;Append(~zv,zz);
 end for;
 for n in[0..4]do for r in[1..Nrows(M0)]do assert &+[&+[RM[j+1][r,c]*zv[n-j+1][c]:c in[1..121]]:j in[0..n]]eq 0;end for;end for;
 return <[Code(Coefficient(xx,n)):n in[0..4]],[Code(Coefficient(yy,n)):n in[0..4]],[[Code(zv[n+1][i]):n in[0..4]]:i in[1..121]]>;
end function;
function RecordPoint(E,P)
 E,P:=SimplifyPoint(E,P);F:=BaseRing(E);z,xy:=CanonicalPoint(E,P);jt:=CanonicalJets(E,P,z);d:=LCM([Degree(MinimalPolynomial(x,GF(2))):x in xy]);assert LCM([Degree(MinimalPolynomial(x,GF(2))):x in z])eq d;
 rec:=<[Integers()!c:c in Eltseq(DefiningPolynomial(F,GF(2)))],[Code(a):a in aInvariants(E)],[Code(P[i]):i in[1,2]],[Code(a):a in z],[Code(a):a in xy],jt>;return rec,d;
end function;

root:="work/deg10/continuation_20261030c/C65/";dat:=eval Read(root cat "source65.m");d,ky:=Explode(dat[3][CASEID]);assert d eq 10;pol:=[r[2]:r in dat[1]|r[1] eq d][1];PF:=PolynomialRing(GF(2));F<a>:=ext<GF(2)|PF!pol>;b:=Evaluate(PF!ky[1..d],a);c:=Evaluate(PF!ky[d+1..2*d],a);ET:=EllipticCurve([F|1-c,-b,-b,0,0]);PT:=ET![0,0];assert Order(PT) eq 65;
ai:=aInvariants(ET);uu:=ai[1];assert uu ne 0;rr:=ai[3]/uu;tr:=(ai[4]+rr^2)/uu;new2:=(ai[2]+rr)/uu^2;new6:=(ai[5]+rr*ai[4]+rr^2*ai[2]+rr^3+tr^2)/uu^6;
E:=EllipticCurve([F|1,new2,0,0,new6]);PP:=E![(PT[1]-rr)/uu^2,(PT[2]-tr)/uu^3];assert Order(PP) eq 65 and TateNormalForm_xy(E,PP) eq TateNormalForm_xy(ET,PT);mp:=MinimalPolynomial(jInvariant(E),GF(2));assert Degree(mp) in [5,10];
SourceRecords:=[* *];SourceHeckeRecords:=[* *];SourceJPolynomials:=[* *];
sr,sd:=RecordPoint(E,PP);assert sd eq 10;Append(~SourceRecords,sr);Append(~SourceJPolynomials,[Integers()!c:c in Eltseq(mp)]);
fast:=MDIsogeniesFast(E,PP,3);assert &+[a[3]:a in fast] eq 4;trs:=[* *];
for dat in fast do
 ET,PT,e:=Explode(dat);assert Order(PT) eq 65;tr,td:=RecordPoint(ET,PT);assert (10*e) mod td eq 0;mult:=10*e div td;Append(~trs,<e,td,mult,tr>);print "C65_DECIC_HECKE_TARGET",CASEID,td,mult,Cputime(tt);
end for;
assert &+[r[2]*r[3]:r in trs] eq 40;Append(~SourceHeckeRecords,trs);
out:=root cat Sprintf("ordinary10_%o.m",CASEID);PrintFile(out,Sprintf("SourceJPolynomials:=%o;SourceRecords:=%o;SourceHeckeRecords:=%o;",SourceJPolynomials,SourceRecords,SourceHeckeRecords):Overwrite:=true);print "C65_DECIC_ELLIPTIC_JETS_COMPLETE",CASEID,sd,[r[2]:r in trs],[r[3]:r in trs],Cputime(tt);quit;
