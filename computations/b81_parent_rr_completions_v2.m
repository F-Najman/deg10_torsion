SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(5*10^9);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");root:="work/deg10/continuation_20261030b/C81/";dat:=eval Read(root cat "y46_kummer_model.m");src:=eval Read(root cat "source81.m");k<eta>:=GF(4);KP:=PolynomialRing(k);B<t>:=FunctionField(k);RF<Y>:=PolynomialRing(B);
function KC(v)return &+[k!v[i]*eta^(i-1):i in[1..#v]];end function;
function Rat(v)return B!(KP![KC(c):c in v[1]])/B!(KP![KC(c):c in v[2]]);end function;
F<y>:=FunctionField(RF![Rat(v):v in dat[1]]);assert Genus(F) eq 13 and Degree(F) eq 8;
function FF(v)return F![Rat(c):c in v];end function;
OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);coords:=[FF(v):v in dat[2]];fk:=FF(dat[3]);nxco:=[FF(v):v in dat[4]];power:=dat[5];AA:=FF(dat[6]);BB:=FF(dat[7]);ac:=FF(dat[8]);Eco:=[FF(v):v in dat[10]];Euniv:=EllipticCurve(Eco);P27:=Euniv![0,0];assert 27*P27 eq Euniv!0 and 9*P27 ne Euniv!0;
T3:=9*P27;a1,a2,a3,a4,a6:=Explode(Eco);xx,yy:=Explode([T3[1],T3[2]]);sl:=(3*xx^2+2*a2*xx+a4-a1*yy)/(2*yy+a1*xx+a3);assert fk eq (P27[2]-yy)-sl*(P27[1]-xx);
RZ<Z>:=PolynomialRing(F);Q:=quo<RZ|Z^3-fk>;xxr:=Q!ac+(Q!AA)*Q.1^power+(Q!BB)*Q.1^(2*power);assert Evaluate(RZ!nxco,xxr) eq Q!0;
function PD(row)
 OO:=row[2] eq 1 select OF else OI;I:=ideal<OO|[FF(v):v in row[3]]>;assert IsPrime(I);p:=Place(I);assert Degree(p) eq row[1];return p;
end function;
ps:=[];vs:=[];for r in dat[9] do p:=PD([*r[1],r[3],r[4]*]);assert Valuation(fk,p) eq r[2];Append(~ps,p);Append(~vs,r[2]);end for;assert Divisor(fk) eq &+[vs[i]*Divisor(ps[i]):i in[1..#ps]];
ram:=[ps[i]:i in[1..#ps]|vs[i] mod 3 ne 0];assert #ram eq 9 and &and[Degree(p) eq 1:p in ram];assert 3*(Genus(F)-1)+1+#ram eq 46;
function Sig(b,c)return [b/(16*b^2+(1-20*c-8*c^2)*b+c*(c-1)^3),b];end function;
sg:=Sig(-Eco[2],1-Eco[1]);b2,c2:=Explode(TateNormalForm_bc(Euniv,2*P27));sg2:=Sig(b2,c2);sigs:=[[Valuation(v,p):v in sg]:p in ram];assert #Seqset(sigs) eq 9;delta:=[Index(sigs,[Valuation(v,p):v in sg2]):p in ram];assert Sort(delta) eq [1..9] and &and[delta[i] ne i:i in[1..9]];
Dzero:=[&+[Floor(Rationals()!(j*vs[i])/3)*Divisor(ps[i]):i in[1..#ps]]:j in[0..2]];
function RankFunctions(bs)
 den:=LCM(&cat[[Denominator(v):v in Eltseq(g)]:g in bs]);cs:=[[KP!(B!den*v):v in Eltseq(g)]:g in bs];dg:=Maximum(&cat[[Degree(v):v in row]:row in cs]);return Rank(Matrix(k,[&cat[[Coefficient(v,i):i in[0..dg]]:v in row]:row in cs]));
end function;
function MapRes(v,info)
 K:=info[4];gv:=info[5];M:=info[6];co:=Eltseq(Solution(M,Vector(k,Eltseq(v,k))));return &+[K!co[i]*gv^(i-1):i in[1..#co]];
end function;
function EV(v,info)
 if v eq 0 then return info[4]!0;end if;assert Valuation(v,info[1]) ge 0;return MapRes(Coefficient(info[3](v),0),info);
end function;
function NormFinite(E,U)
 ai:=aInvariants(E);T:=3*U;x:=T[1];yy:=T[2];a3:=2*yy+ai[1]*x+ai[3];a2:=3*x+ai[2];a4:=3*x^2+2*ai[2]*x+ai[4]-ai[1]*yy;sl:=a4/a3;a2-:=sl^2+ai[1]*sl;kap:=a3/a2;b:=-a2^3/a3^2;c:=1-(ai[1]+2*sl)*a2/a3;ET:=EllipticCurve([BaseRing(E)|1-c,-b,-b,0,0]);V:=ET![(U[1]-x)/kap^2,(U[2]-yy-sl*(U[1]-x))/kap^3];assert 3*V eq ET![0,0];return ET,V;
end function;
function LocalInfo(lb,ci)
 p:=PD(lb[1]);mult:=lb[2];tau:=FF(lb[3]);assert Valuation(tau,p) eq 1;K<z>:=ext<k|KP![KC(c):c in lb[4][1]]>;assert Degree(K,k) eq Degree(p);d,ky:=Explode(src[3][ci]);b:=Evaluate(KP!ky[1..d],z);c:=Evaluate(KP!ky[d+1..2*d],z);E:=EllipticCurve([K|1-c,-b,-b,0,0]);P:=E![0,0];assert 81*P eq E!0 and 27*P ne E!0;U:=mult eq 1 select P else 2*P;vals:=TateNormalForm_xy(E,3*U);
 S,ex:=Completion(F,p:Precision:=24);ts:=ex(tau);assert Valuation(ts) eq 1 and AbsolutePrecision(ts) ge 4;u1:=Coefficient(ts,1);assert u1 ne 0;found:=false;
 for bits in [[k|1,0],[k|0,1]] cat [[k|1,a]:a in k|a ne 0] do
  gf:=&+[F!bits[i]*coords[i]:i in[1,2]];gc:=Coefficient(ex(gf),0);po:=MinimalPolynomial(gc,k);if Degree(po) ne Degree(p) or Evaluate(po,gf) ne tau then continue;end if;gv:=&+[K!bits[i]*vals[i]:i in[1,2]];if Evaluate(po,gv) ne 0 then continue;end if;
  M:=Matrix(k,[Eltseq(gc^i,k):i in[0..Degree(p)-1]]);assert Rank(M) eq Degree(p);info:=[*p,S,ex,K,gv,M,tau,u1*];if [EV(v,info):v in coords] ne vals then continue;end if;found:=true;break;
 end for;assert found;
 function DK(v)return &+[K!KC(v[i])*z^(i-1):i in[1..#v]];end function;
 roots:=[DK(v):v in lb[5]];actual:=DK(lb[6]);assert #Seqset(roots) eq 3 and actual in roots and &and[w^3 eq EV(fk,info):w in roots];ET,V:=NormFinite(E,U);assert aInvariants(ET) eq [EV(v,info):v in Eco];S3:=9*(ET![0,0]);RX<X>:=PolynomialRing(K);EG,phi:=IsogenyFromKernel(ET,X-S3[1]);r:=phi(V)[1];assert Evaluate(RX![EV(v,info):v in nxco],r) eq 0;matches:=[w:w in roots|EV(ac,info)+EV(AA,info)*w^power+EV(BB,info)*w^(2*power) eq r];assert matches eq [actual];
 return info,roots,actual;
end function;
infos:=AssociativeArray();rootscache:=AssociativeArray();actualcache:=AssociativeArray();result:=[* *];nspaces:=0;nlocal:=0;
for ci in[1,2,3,7,10] do
 d:=src[3][ci][1];csidx:=d eq 9 select [1..9] else [0];
 for cidx in csidx do
  p:=eval Read(root cat Sprintf("rr_case%o_cusp%o.m",ci,cidx));NC:=DivisorGroup(F)!0;
  for lb in p[3] do ky:=Sprint([*ci,lb*]);if not IsDefined(infos,ky) then info,rr,act:=LocalInfo(lb,ci);infos[ky]:=info;rootscache[ky]:=rr;actualcache[ky]:=act;nlocal+:=1;end if;NC+:=lb[2]*Divisor(infos[ky][1]);end for;
  assert Degree(NC) eq 3*d;divs:=[NC+z:z in Dzero];
  if cidx gt 0 then for j in[0..2] do divs[j+1]+:=(Floor(Rationals()!(1+j*Valuation(fk,ram[cidx]))/3)-Floor(Rationals()!(j*Valuation(fk,ram[cidx]))/3))*Divisor(ram[cidx]);divs[j+1]+:=(Floor(Rationals()!(2+j*Valuation(fk,ram[delta[cidx]]))/3)-Floor(Rationals()!(j*Valuation(fk,ram[delta[cidx]]))/3))*Divisor(ram[delta[cidx]]);end for;end if;
  entries:=[* *];dims:=[];
  for jj in[1..3] do assert p[1][jj][1] eq jj-1;bs:=[FF(v):v in p[1][jj][2]];db:=#Basis(divs[jj]:Reduction:=false,Simplification:="None");assert #bs eq db and RankFunctions(bs) eq db;I,J:=Ideals(divs[jj]);assert &and[IsIntegral(I*(g*OF)) and IsIntegral(J*(g*OI)):g in bs];Append(~dims,db);nspaces+:=1;for g in bs do Append(~entries,[*jj-1,g*]);end for;end for;
  rows:=[[k|]:i in[1..#entries]];
  for lb in p[3] do ky:=Sprint([*ci,lb*]);info:=infos[ky];mult:=lb[2];tau:=info[7];rr:=rootscache[ky];actual:=actualcache[ky];others:=[w:w in rr|w ne actual];fp:=MapRes(Coefficient(info[3](fk),1)/info[8],info);jets:=[* *];
   for entry in entries do g:=tau^mult*entry[2];assert Valuation(g,info[1]) ge 0;gs:=info[3](g);assert AbsolutePrecision(gs) ge 2;Append(~jets,[*MapRes(Coefficient(gs,0),info),MapRes(Coefficient(gs,1)/info[8],info)*]);end for;
   for w in others do for ord in[0..mult-1] do for i in[1..#entries] do j:=entries[i][1];v:=ord eq 0 select w^j*jets[i][1] else w^j*jets[i][2]+j*w^(j-1)*fp/w^2*jets[i][1];rows[i] cat:=Eltseq(v,k);end for;end for;end for;
  end for;
  MM:=Matrix(k,rows);old:=Matrix(k,[[KC(v):v in row]:row in p[2]]);assert MM eq old;rank:=Rank(MM);assert #entries-rank eq 1 and Ncols(MM) eq 6*d;assert dims eq (d eq 9 select[15,11,13]else[18,12,15]);Append(~result,[*ci,cidx,dims,Nrows(MM),Ncols(MM),rank,1*]);print "C81_PARENT_COMPLETION_RR_CASE",ci,cidx,dims,rank,Cputime(tt);
 end for;
end for;
assert #result eq 29 and nspaces eq 87 and nlocal eq 14;Write(root cat "parent_native_completion_RR_v2.m",Sprint([*result,nspaces,nlocal*]):Overwrite:=true);print "C81_PARENT_COMPLETION_RR_V2_COMPLETE",29,nspaces,nlocal,Cputime(tt);quit;
