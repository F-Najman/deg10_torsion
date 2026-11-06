SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(5*10^9);tt:=Cputime();AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");root:="work/deg10/continuation_20261030b/C81/";dat:=eval Read(root cat "y46_kummer_model.m");k<eta>:=GF(4);KP:=PolynomialRing(k);B<t>:=FunctionField(k);RF<Y>:=PolynomialRing(B);
function KC(v) return &+[k!v[i]*eta^(i-1):i in [1..#v]];end function;
function Rat(v) return B!(KP![KC(c):c in v[1]])/B!(KP![KC(c):c in v[2]]);end function;
f:=RF![Rat(v):v in dat[1]];f/:=LeadingCoefficient(f);F<y>:=FunctionField(f);assert Genus(F) eq 13;
function FF(v) return F![Rat(w):w in v];end function;
coords:=[FF(v):v in dat[2]];fk:=FF(dat[3]);nxco:=[FF(v):v in dat[4]];power:=dat[5];AA:=FF(dat[6]);BB:=FF(dat[7]);a2c:=FF(dat[8]);Eco:=[FF(v):v in dat[10]];Euniv:=EllipticCurve(Eco);P27:=Euniv![0,0];OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);ps:=[];vs:=[];
for row in dat[9] do OO:=row[3] eq 1 select OF else OI;I:=ideal<OO|[FF(v):v in row[4]]>;p:=Place(I);assert Degree(p) eq row[1] and Valuation(fk,p) eq row[2];Append(~ps,p);Append(~vs,row[2]);end for;assert Divisor(fk) eq &+[vs[i]*Divisor(ps[i]):i in [1..#ps]];
ram:=[p:p in ps|Valuation(fk,p) mod 3 ne 0];assert #ram eq 9 and &and[Degree(p) eq 1:p in ram];
function Sig(b,c) return [b/(16*b^2+(1-20*c-8*c^2)*b+c*(c-1)^3),b];end function;
b:=-Eco[2];c:=1-Eco[1];sg:=Sig(b,c);b2,c2:=Explode(TateNormalForm_bc(Euniv,2*P27));sg2:=Sig(b2,c2);signatures:=[[Valuation(v,p):v in sg]:p in ram];assert #Seqset(signatures) eq 9;delta:=[];for p in ram do imsig:=[Valuation(v,p):v in sg2];j:=Index(signatures,imsig);assert j gt 0;Append(~delta,j);end for;assert Sort(delta) eq [1..9] and &and[delta[i] ne i:i in [1..9]];
// Exact derivation with respect to the rational base parameter.
function DR(v) n:=Numerator(v);d:=Denominator(v);return B!(Derivative(n)*d-n*Derivative(d))/B!(d^2);end function;
yprime:=-Evaluate(RF![DR(v):v in Coefficients(f)],F.1)/Evaluate(Derivative(f),F.1);
function D(v) cs:=Eltseq(F!v);return &+[F!DR(cs[i])*F.1^(i-1)+(i eq 1 select F!0 else (i-1)*(F!cs[i])*F.1^(i-2)*yprime):i in [1..#cs]];end function;
assert D(F!B.1) eq 1;
function PlaceFromPair(vals)
 K:=Parent(vals[1]);e:=LCM([Degree(MinimalPolynomial(v,k)):v in vals]);found:=false;
 for bits in [[k|1,0],[k|0,1]] cat [[k|1,a]:a in k|a ne 0] do gen:=&+[K!bits[i]*vals[i]:i in [1,2]];po:=MinimalPolynomial(gen,k);if Degree(po) eq e then chosen:=bits;found:=true;break;end if;end for;assert found;
 MM:=Matrix(k,[Eltseq(gen^i,k):i in [0..e-1]]);assert Rank(MM) eq e;cc:=[Eltseq(Solution(MM,Vector(k,Eltseq(v,k)))):v in vals];gg:=&+[F!chosen[i]*coords[i]:i in [1,2]];fs:=[Evaluate(po,gg)] cat [coords[j]-&+[F!cc[j][i+1]*gg^i:i in [0..e-1]]:j in [1,2]];ans:=[];
 for OO in [OF,OI] do II:=ideal<OO|fs>;for fac in Factorization(II) do if fac[2] le 0 then continue;end if;pl:=Place(fac[1]);if Degree(pl) eq e and &and[g eq 0 or Valuation(g,pl) gt 0:g in fs] and not pl in ans then Append(~ans,pl);end if;end for;end for;assert #ans eq 1;return ans[1];
end function;
function Normalize(E,U)
 ai:=aInvariants(E);T:=3*U;x:=T[1];yy:=T[2];a3:=2*yy+ai[1]*x+ai[3];a2:=3*x+ai[2];a4:=3*x^2+2*ai[2]*x+ai[4]-ai[1]*yy;s:=a4/a3;a2-:=s^2+ai[1]*s;kap:=a3/a2;b:=-a2^3/a3^2;c:=1-(ai[1]+2*s)*a2/a3;ET:=EllipticCurve([BaseRing(E)|1-c,-b,-b,0,0]);V:=ET![(U[1]-x)/kap^2,(U[2]-yy-s*(U[1]-x))/kap^3];assert 3*V eq ET![0,0];return ET,V;
end function;
function EV(v,loc)
 if v eq 0 then return loc[2]!0;end if;pl:=loc[1];K:=loc[2];gr:=loc[4];kr:=Parent(gr);z:=kr!Evaluate(v,pl);cs:=Eltseq(Solution(loc[5],Vector(k,Eltseq(z,k))));return &+[K!cs[i]*loc[3]^(i-1):i in [1..#cs]];
end function;
function LocalInput(K,E,U,mult)
 vals:=TateNormalForm_xy(E,3*U);pl:=PlaceFromPair(vals);e:=Degree(pl);assert e eq Degree(K,k) and Valuation(fk,pl) eq 0;found:=false;
 for bits in [[k|1,0],[k|0,1]] cat [[k|1,a]:a in k|a ne 0] do
  gv:=&+[K!bits[i]*vals[i]:i in [1,2]];po:=MinimalPolynomial(gv,k);if Degree(po) ne e then continue;end if;gf:=&+[F!bits[i]*coords[i]:i in [1,2]];tau:=Evaluate(po,gf);if Valuation(tau,pl) ne 1 then continue;end if;found:=true;break;
 end for;assert found;gr:=Evaluate(gf,pl);assert MinimalPolynomial(gr,k) eq po;MM:=Matrix(k,[Eltseq(gr^i,k):i in [0..e-1]]);assert Rank(MM) eq e;loc:=[*pl,K,gv,gr,MM,tau,[K|],K!0,mult*];
 ET,V:=Normalize(E,U);assert [EV(v,loc):v in Eco] eq aInvariants(ET);P:=ET![0,0];S:=9*P;RX<X>:=PolynomialRing(K);EG,phi:=IsogenyFromKernel(ET,X-S[1]);r:=phi(V)[1];nc:=[EV(v,loc):v in nxco];assert Evaluate(RX!nc,r) eq 0;fv:=EV(fk,loc);rr:=Roots(X^3-fv);assert #rr eq 3;roots:=[q[1]:q in rr];av:=EV(a2c,loc);Aa:=EV(AA,loc);Bb:=EV(BB,loc);choices:=[w:w in roots|av+Aa*w^power+Bb*w^(2*power) eq r];assert #choices eq 1;loc[7]:=roots;loc[8]:=choices[1];return loc;
end function;
function EC(v) return [Integers()!a:a in Eltseq(k!v,GF(2))];end function;
function ER(v) return [[EC(c):c in Coefficients(Numerator(v))],[EC(c):c in Coefficients(Denominator(v))]];end function;
function EF(v) return [ER(c):c in Eltseq(F!v)];end function;
function PD(p) return [*Degree(p),Order(Ideal(p)) eq OF select 1 else 0,[EF(v):v in Basis(Ideal(p))]*];end function;
source:=eval Read(root cat "source81.m");d,ky:=Explode(source[3][CASEID]);assert d in [9,10];pol:=[r[2]:r in source[1]|r[1] eq d][1];factors:=Factorization(KP!pol);assert &and[q[2] eq 1:q in factors];localdata:=[* *];
for fac in factors do
 K<z>:=ext<k|fac[1]>;b:=Evaluate(KP!ky[1..d],z);c:=Evaluate(KP!ky[d+1..2*d],z);E:=EllipticCurve([K|1-c,-b,-b,0,0]);P:=E![0,0];assert 81*P eq E!0 and 27*P ne E!0;
 Append(~localdata,LocalInput(K,E,P,1));Append(~localdata,LocalInput(K,E,2*P,2));
end for;
assert #Seqset([v[1]:v in localdata]) eq #localdata;assert &+[Degree(v[1])*v[9]:v in localdata] eq 3*d;
NC:=&+[v[9]*Divisor(v[1]):v in localdata];Dzero:=[&+[Floor(Rationals()!(j*vs[i])/3)*Divisor(ps[i]):i in [1..#ps]]:j in [0,1,2]];print "C81_RR_LOCAL_INPUTS_READY",CASEID,d,#localdata,delta,Cputime(tt);
csidx:=d eq 9 select [1..9] else [0];results:=[* *];
for ci in csidx do
 divs:=[NC+Dzero[j+1]:j in [0,1,2]];
 if ci gt 0 then
  assert delta[ci] ne ci;
  for j in [0,1,2] do divs[j+1]+:=(Floor(Rationals()!(1+j*Valuation(fk,ram[ci]))/3)-Floor(Rationals()!(j*Valuation(fk,ram[ci]))/3))*Divisor(ram[ci]);divs[j+1]+:=(Floor(Rationals()!(2+j*Valuation(fk,ram[delta[ci]]))/3)-Floor(Rationals()!(j*Valuation(fk,ram[delta[ci]]))/3))*Divisor(ram[delta[ci]]);end for;
 end if;
 bases:=[*Basis(v:Reduction:=false,Simplification:="None"):v in divs*];entries:=[* *];for j in [0..2] do for v in bases[j+1] do Append(~entries,[*j,v*]);end for;end for;nr:=#entries;rows:=[ [k|]:i in [1..nr] ];
 for loc in localdata do
  pl:=loc[1];tau:=loc[6];dtau:=D(tau);assert dtau ne 0;mult:=loc[9];K:=loc[2];wp:=loc[8];others:=[w:w in loc[7]|w ne wp];assert #others eq 2;fprime:=EV(D(fk)/dtau,loc);
  gg:=[* *];for entry in entries do j:=entry[1];g:=tau^mult*entry[2];assert Valuation(g,pl) ge 0;dg:=D(g)/dtau;assert Valuation(dg,pl) ge 0;Append(~gg,[*EV(g,loc),EV(dg,loc)*]);end for;
  for w in others do
   for order in [0..mult-1] do
    for i in [1..nr] do j:=entries[i][1];val:=order eq 0 select w^j*gg[i][1] else w^j*gg[i][2]+j*w^(j-1)*fprime/w^2*gg[i][1];rows[i] cat:=Eltseq(val,k);end for;
   end for;
  end for;
 end for;
 MM:=Matrix(k,rows);rank:=Rank(MM);dim:=nr-rank;assert dim ge 1;
 proof:=[*[ [*j,[EF(v):v in bases[j+1]]*]:j in [0,1,2]],[[EC(v):v in row]:row in rows],[[*PD(loc[1]),loc[9],EF(loc[6]),[[EC(v):v in Coefficients(DefiningPolynomial(loc[2],k))]],[[EC(v):v in Eltseq(w,k)]:w in loc[7]],[EC(v):v in Eltseq(loc[8],k)]*]:loc in localdata]*];
 Write(root cat Sprintf("rr_case%o_cusp%o.m",CASEID,ci),Sprint(proof):Overwrite:=true);Append(~results,[*CASEID,ci,[#v:v in bases],nr,Ncols(MM),rank,dim*]);print "C81_Y46_RR_JET_TEST",CASEID,ci,[#v:v in bases],nr,Ncols(MM),rank,dim,Cputime(tt);
end for;
Write(root cat Sprintf("rr_case%o_summary.m",CASEID),Sprint(results):Overwrite:=true);print "C81_Y46_RR_JET_COMPLETE",CASEID,#results,[r[7]:r in results],Cputime(tt);quit;
