SetColumns(0);SetNthreads(1);SetQuitOnError(true);SetMemoryLimit(8*10^9);Alarm(800);tt:=Cputime();root:="work/deg10/continuation_20261031b/M3N33/";dat:=eval Read(root cat "relative_kummer_frame.m");bps:=eval Read(root cat "relative_base_cusps.m");F2:=GF(2);P2:=PolynomialRing(F2);k<w>:=ext<F2|P2!dat[1]>;B<u>:=FunctionField(k);PK:=PolynomialRing(k);RY<Y>:=PolynomialRing(B);
function DC(a)return k!(a mod 2)+k!(a div 2)*w;end function;
function EC(a)v:=Eltseq(k!a,F2);return Integers()!v[1]+2*Integers()!v[2];end function;
function DR(a)return B!(PK![DC(c):c in a[1]])/(B!(PK![DC(c):c in a[2]]));end function;
F<v>:=FunctionField(RY![DR(c):c in dat[2]]);assert Genus(F)eq 21;
function DF(a)return F![DR(c):c in a];end function;
OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
function BasePrime(a)O:=a[2]eq 1 select OF else OI;I:=ideal<O|[DF(c):c in a[3]]>;assert IsPrime(I);P:=Place(I);assert Degree(P)eq a[1];return P;end function;
dd:=-DF(dat[3][1]);assert DF(dat[3][2])eq 0 and DF(dat[3][3])eq 0 and DF(dat[3][4])eq 1;TM:=Matrix(k,81,81,[DC(c):c in dat[8][2]]);U:=Matrix(k,#dat[8][3]div 81,81,[DC(c):c in dat[8][3]]);assert Nrows(U)eq 36 and Rank(U)eq 36;
functions:=[F!u,v,dd];idx:=[* *];for bf in dat[5]do ids:=[];for co in bf do f:=DF(co);i:=Index(functions,f);if i eq 0 then Append(~functions,f);i:=#functions;end if;Append(~ids,i);end for;Append(~idx,ids);end for;
function ES(z,K)
 if z eq 0 then return[*0,140,[]*];end if;a:=Valuation(z);b:=AbsolutePrecision(z);return[*a,b,[[EC(c):c in Eltseq(K!Coefficient(z,j),k)]:j in[a..b-1]]*];
end function;
function SubstituteMonomial(A,e,c)
 C:=Parent(A);s:=C.1;if A eq 0 then return C!0;end if;lo:=Valuation(A);hi:=AbsolutePrecision(A);ans:=&+[C|Coefficient(A,j)*c^j*s^(e*j):j in[lo..hi-1]];return ans+O(s^(e*hi));
end function;
function CubeRootUnit(V,beta)
 C:=Parent(V);s:=C.1;assert Valuation(V)eq 0 and beta^3 eq Coefficient(V,0);N:=AbsolutePrecision(V);z:=C!beta+O(s^N);for j in[1..12]do z:=V/z^2;if z^3 eq V then break;end if;end for;assert z^3 eq V;return z;
end function;
rows:=[* *];wits:=[* *];for ii in[1..#bps]do
 P:=BasePrime(bps[ii]);d:=Degree(P);assert d in[1,5];CP,cp:=Completion(F,P:Precision:=140);K:=BaseRing(CP);Embed(k,K);s:=CP.1;series:=[cp(f):f in functions];assert &and[AbsolutePrecision(a)ge 100:a in series];ds:=series[3];m:=Valuation(ds);V:=ds/s^m;v0:=Coefficient(V,0);assert v0 ne 0;
 if m mod 3 eq 0 then e:=1;c:=K!1;PP:=PolynomialRing(K);betas:=[r[1]:r in Roots(PP.1^3-v0)];assert #betas eq 3;else e:=3;a:=m mod 3;b:=(1-a*m)div 3;c:=v0^(-a);betas:=[v0^b];assert c^m*v0 eq betas[1]^3;end if;
 for jj in[1..#betas]do
  beta:=betas[jj];ss:=[SubstituteMonomial(z,e,c):z in series];VV:=c^m*SubstituteMonomial(V,e,c);H:=CubeRootUnit(VV,beta);rr:=s^(e*m div 3)*H;assert rr^3 eq ss[3];du:=Derivative(ss[1]);assert du ne 0;hol:=[(&+[ss[ids[j]]*rr^(j-1):j in[1..#ids]])*du:ids in idx];assert &and[z eq 0 or Valuation(z)ge 0:z in hol];assert &and[AbsolutePrecision(z)ge 5:z in hol];nv:=d eq 1 select 5 else 1;J:=Matrix(k,81,nv*d,&cat[&cat[Eltseq(K!Coefficient(z,j),k):j in[0..nv-1]]:z in hol]);assert Rank(J)eq nv*d;UJ:=U*J;
  Append(~rows,[*ii,jj,d,e,m,[EC(a):a in Eltseq(beta,k)],[EC(a):a in Eltseq(c,k)],[[EC(a):a in row]:row in Rows(J)],[[EC(a):a in row]:row in Rows(UJ)]*]);print "M333_ANALYTIC_CUSP_JET",ii,jj,d,e,m,Rank(UJ),Cputime(tt);
 end for;
 Append(~wits,[*ii,d,m,[EC(a):a in Coefficients(DefiningPolynomial(K,k))],[ES(z,K):z in series]*]);Write(root cat "analytic_cusp_jets_v2_progress.m",Sprint([*rows,wits*]):Overwrite:=true);
end for;
assert #rows eq 48 and #[r:r in rows|r[3]eq 1]eq 40 and #[r:r in rows|r[3]eq 5]eq 8;perm:=[];spaces:=[RowSpace(Transpose(Matrix(k,81,r[3],&cat[[DC(c):c in row[1..r[3]]]:row in r[8]]))):r in rows];
for i in[1..#rows]do J:=Matrix(k,81,rows[i][3],&cat[[DC(c):c in row[1..rows[i][3]]]:row in rows[i][8]]);sp:=RowSpace(Transpose(TM*J));hits:=[j:j in[1..#rows]|rows[j][3]eq rows[i][3]and spaces[j]eq sp];assert #hits eq 1 and hits[1]ne i;Append(~perm,hits[1]);end for;assert &and[perm[perm[i]]eq i:i in[1..#rows]];
Write(root cat "analytic_cusp_jets_v2.m",Sprint([*rows,wits,perm,idx*]):Overwrite:=true);print "D31B_M333_ANALYTIC_CUSP_JETS_COMPLETE",#rows,perm,Cputime(tt);quit;
