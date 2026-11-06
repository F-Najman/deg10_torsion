SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);SetNthreads(1);Alarm(600);tt:=Cputime();root:="work/deg10/continuation_20261031b/C66/";
load "work/deg10/continuation_20261031b/C66/inputs/F5_rr/input33.m";
k:=GF(5);K<t>:=FunctionField(k);PK:=PolynomialRing(k);RY<y>:=PolynomialRing(K);
function Rat(v)return (K!PK!v[1])/(K!PK!v[2]);end function;
F<u>:=FunctionField(RY![Rat(v):v in FieldEquation]);assert Degree(F)eq 12 and Genus(F)eq 21;
function Fun(v)return F![Rat(w):w in v];end function;
OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);blocks:=[* *];
for z in GlobalCuspBlocks do
 sig,deg,pp:=Explode(z);DD:=DivisorGroup(F)!0;
 for v in pp do dd,m,flag,ib:=Explode(v);O:=flag eq 1 select OF else OI;II:=ideal<O|[Fun(b):b in ib]>;assert IsPrime(II);P:=Place(II);assert Degree(P)eq dd;DD+:=m*Divisor(P);end for;
 assert Degree(DD)eq deg;Append(~blocks,DD);
end for;

load "work/deg10/continuation_20261031b/C66/gonal_cusp_input.m";
function DPole(f)return &+([DivisorGroup(F)!0] cat [-Valuation(f,p)*Divisor(p):p in Poles(f)]);end function;
function DZero(f)return &+([DivisorGroup(F)!0] cat [Valuation(f,p)*Divisor(p):p in Zeros(f)]);end function;
x:=F.1;y:=F!K.1;r:=(x^2*y-x*y+y-1)/x/(x*y-1);ss0:=(x*y-y+1)/x/y;b:=r*ss0*(r-1);c:=ss0*(r-1);E:=EllipticCurve([1-c,-b,-b,F!0,F!0]);P:=E![F|0,0];assert 33*P eq E!0 and 3*P ne E!0 and 11*P ne E!0;units:=[a:a in[1..16]|GCD(a,33) eq 1];xy:=[x,y];
D0:=&+[GonalCuspVectors[1][i]*blocks[i]:i in[1..17]];VV,vm:=RiemannRochSpace(D0);assert Dimension(VV) eq 2;fs:=[vm(w):w in Basis(VV)];g0:=[f:f in fs|Degree(DPole(f)) gt 0][1];assert Degree(DPole(g0)) eq 10;
gs:=[F|];for a in units do xy0:=TateNormalForm_xy(E,a*P);ss:=hom<BaseRing(F)->F|xy0[2]>;dd:=hom<F->F|ss,xy0[1]>;Append(~gs,dd(g0));end for;assert &and[Degree(DPole(g)) eq 10:g in gs];j:=jInvariant(E);cs:=Poles(j);assert &and[Seqset(Poles(g)) subset Seqset(cs):g in gs];
function ER(f)return [[Integers()!z:z in Coefficients(Numerator(f))],[Integers()!z:z in Coefficients(Denominator(f))]];end function;
function EF(f)return [ER(z):z in Eltseq(F!f)];end function;
function PD(p)g,h:=TwoGenerators(p);OO:=IsFinite(p) select MaximalOrderFinite(F) else MaximalOrderInfinite(F);assert ideal<OO|g,h> eq Ideal(p);return[*Degree(p),IsFinite(p) select 1 else 0,[EF(g),EF(h)]*];end function;
function CV(w)return [Integers()!q:q in Eltseq(w,GF(5))];end function;
function Enc(E,P)return[*[CV(a):a in aInvariants(E)],[CV(P[1]),CV(P[2])]*];end function;
Write(root cat "gonal_functions_F5.m",Sprint([*units,[EF(g):g in gs],[EF(a):a in aInvariants(E)],xy[1] eq F!K.1,[PD(p):p in cs]*]):Overwrite:=true);
vectors:=[];for g in gs do DD:=DPole(g);vv:=[Valuation(DD,Support(B)[1]):B in blocks];assert DD eq &+[vv[i]*blocks[i]:i in[1..17]];Append(~vectors,vv);end for;Write(root cat "pencil_cusp_vectors_F5.m",Sprint([*[b[1]:b in GlobalCuspBlocks],units,vectors*]):Overwrite:=true);printf "D31B_C66_GONAL_MODEL_READY functions10degree10=true time=%o\n",Cputime(tt);out:=[* *];g:=gs[1];assert units[1] eq 1;allcounts:=[];
for aa in[0..5] do
 D:=aa eq 5 select DPole(g) else DZero(g-GF(5)!aa);assert Degree(D) eq 10;ps,es:=Support(D);bounds:=[0:i in[1..60]];records:=[* *];cuspmass:=0;
 for ii in[1..#ps] do
  p:=ps[ii];m:=es[ii];d:=Degree(p);
  if p in cs then cuspmass+:=m*d;continue;end if;
  ai:=[Evaluate(v,p):v in aInvariants(E)];L:=Universe(ai);EE:=EllipticCurve(ai);PP:=EE![L|0,0];assert Order(PP) eq 33;RX<T>:=PolynomialRing(L);a1,a2,a3,a4,a6:=Explode(ai);cub:=4*T^3+(a1^2+4*a2)*T^2+2*(a1*a3+2*a4)*T+(a3^2+4*a6);assert Degree(cub) eq 3 and IsSquarefree(cub);facs:=Factorization(cub);assert &+[Degree(h[1])*h[2]:h in facs] eq 3;
  for ff in facs do
   assert ff[2] eq 1;e:=Degree(ff[1]);if e eq 1 then M:=L;r:=Roots(ff[1])[1][1];else M<r>:=ext<L|ff[1]>;end if;
   EG:=ChangeRing(EE,M);PG:=EG!PP;R:=EG![r,-(M!a1*r+M!a3)/2];assert Order(R) eq 2;E2,phi:=IsogenyFromKernel(EG,PolynomialRing(M)![-r,1]);Q:=phi(PG);assert Order(Q) eq 33;vals:=[* *];for a in units do xy0:=TateNormalForm_xy(E2,a*Q);zz0:=xy0[1];val:=&+[Evaluate(Eltseq(g0)[i],xy0[2])*zz0^(i-1):i in[1..#Eltseq(g0)]];Append(~vals,val);end for;bs:=[];
   for k in[1..10] do
    b:=vals[k];if b^5 eq b then bb:=Integers()!(GF(5)!b);bounds[(k-1)*6+bb+1]+:=m*d*e;Append(~bs,[k,bb]);end if;
   end for;
   pol:=[Integers()!a:a in Coefficients(DefiningPolynomial(M,GF(5)))];Append(~records,[*PD(p),m,d,e,pol,Enc(EG,PG),CV(r),Enc(E2,Q),bs*]);
  end for;
 end for;
 Append(~out,[*aa,cuspmass,bounds,records*]);Append(~allcounts,Max(bounds));Write(root cat "gonal_fibre_filter_progress.m",Sprint(out):Overwrite:=true);printf "D31B_C66_GONAL_FIBRE a=%o cuspmass=%o maxNCpairmass=%o time=%o\n",aa,cuspmass,Max(bounds),Cputime(tt);
end for;
Write(root cat "gonal_fibre_filter.m",Sprint(out):Overwrite:=true);printf "D31B_C66_GONAL_FIBRE_FILTER_COMPLETE maxima=%o time=%o\n",allcounts,Cputime(tt);quit;
