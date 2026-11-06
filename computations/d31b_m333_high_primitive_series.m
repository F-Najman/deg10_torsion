SetColumns(0);SetNthreads(1);SetQuitOnError(true);SetMemoryLimit(8*10^9);Alarm(800);tt:=Cputime();root:="work/deg10/continuation_20261031b/M3N33/";dat:=eval Read(root cat "relative_kummer_frame.m");bps:=eval Read(root cat "relative_base_cusps.m");low:=eval Read(root cat "analytic_cusp_jets_v4.m");F2:=GF(2);P2:=PolynomialRing(F2);k<w>:=ext<F2|P2!dat[1]>;B<u>:=FunctionField(k);PK:=PolynomialRing(k);RY<Y>:=PolynomialRing(B);
function DC(a)return k!(a mod 2)+k!(a div 2)*w;end function;
function EC(a)v:=Eltseq(k!a,F2);return Integers()!v[1]+2*Integers()!v[2];end function;
function DR(a)return B!(PK![DC(c):c in a[1]])/(B!(PK![DC(c):c in a[2]]));end function;
F<v>:=FunctionField(RY![DR(c):c in dat[2]]);assert Genus(F)eq 21;
function DF(a)return F![DR(c):c in a];end function;
OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
function BasePrime(a)O:=a[2]eq 1 select OF else OI;I:=ideal<O|[DF(c):c in a[3]]>;assert IsPrime(I);P:=Place(I);assert Degree(P)eq a[1];return P;end function;
function ES(z,K)
 if z eq 0 then return[*0,1000,[]*];end if;a:=Valuation(z);b:=AbsolutePrecision(z);return[*a,b,[[EC(c):c in Eltseq(K!Coefficient(z,j),k)]:j in[a..b-1]]*];
end function;
out:=[* *];for i in[1..#bps]do
 P:=BasePrime(bps[i]);old:=low[2][i];assert old[1]eq i and old[2]eq Degree(P);CP,cp:=Completion(F,P:Precision:=1000);K:=BaseRing(CP);Embed(k,K);kpol:=[EC(a):a in Coefficients(DefiningPolynomial(K,k))];assert kpol eq old[4];series:=[cp(F!u),cp(v)];assert &and[AbsolutePrecision(z)ge 900:z in series];
 for j in[1,2]do
  z:=series[j];a,b,cs:=Explode(old[5][j]);if #cs eq 0 then assert z eq 0 or Valuation(z)ge b;else assert #cs eq b-a and Valuation(z)ge a;
   for n in[a..b-1]do vals:=cs[n-a+1];c:=&+[K|DC(vals[m])*K.1^(m-1):m in[1..#vals]];assert Coefficient(z,n)eq c;end for;
  end if;
 end for;
 Append(~out,[*i,Degree(P),kpol,[ES(z,K):z in series]*]);Write(root cat "high_primitive_series_progress.m",Sprint(out):Overwrite:=true);print "M333_HIGH_PRIMITIVES_BASE",i,Degree(P),Cputime(tt);
end for;assert #out eq 24;Write(root cat "high_primitive_series.m",Sprint(out):Overwrite:=true);print "D31B_M333_HIGH_PRIMITIVES_COMPLETE",#out,Cputime(tt);quit;
