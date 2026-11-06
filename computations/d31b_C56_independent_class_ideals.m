SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(6*10^9);Alarm(550);tt:=Cputime();AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");root:="work/deg10/continuation_20261031b/C56/";a:=eval Read(root cat "inputs/target28_proven_picard.m");X:=MDX1(28,GF(3));C:=Curve(X);F0,mp:=AlgorithmicFunctionField(FunctionField(C));F:=RationalExtensionRepresentation(F0);assert Genus(F)eq 10;k:=GF(3);PK:=PolynomialRing(k);B:=BaseField(F);
function Rat(v)return B!(PK!v[1])/B!(PK!v[2]);end function;
assert Coefficients(DefiningPolynomial(F))eq[Rat(v):v in a[1]];
function FF(v)return F![Rat(z):z in v];end function;
xy:=[F!mp(v):v in X`_coordinates_xy];assert xy eq[FF(v):v in a[7]];OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
function DecodePrime(v)
 O:=v[2]eq 1 select OF else OI;J:=ideal<O|[FF(z):z in v[3]]>;assert IsPrime(J);P:=Place(J);assert Degree(P)eq v[1];return P;
end function;
fb:=[DecodePrime([*a[3][i],a[2][i],a[4][i]*]):i in[1..#a[3]]];AB:=DecodePrime(a[12]);assert Degree(AB)eq 1;A:=Divisor(AB);ia:=Index(fb,AB);assert ia gt 0;Z:=FreeAbelianGroup(#a[6][1]);G,q:=quo<Z|[Z!v:v in a[6]]>;fc:=[q(Z!v):v in a[5]];
packet:=eval Read(root cat "inputs/fourmap_filter.m");fields:=AssociativeArray();for r in packet[1]do fields[r[1]]:=ext<k|PK!r[2]>;end for;
function FromKey(key)
 d:=key[1];K:=fields[d];v:=[Evaluate(PK!key[2..d+1],K.1),Evaluate(PK!key[d+2..2*d+1],K.1)];e:=LCM([Degree(MinimalPolynomial(z,k)):z in v]);found:=false;
 for bits in[[1,0],[0,1],[1,1],[1,2]]do gen:=&+[K!bits[j]*v[j]:j in[1,2]];pol:=MinimalPolynomial(gen,k);if Degree(pol)eq e then found:=true;chosen:=bits;break;end if;end for;assert found;
 M:=Matrix(k,[Eltseq(gen^i,k):i in[0..e-1]]);assert Rank(M)eq e;cs:=[Eltseq(Solution(M,Vector(k,Eltseq(z,k)))):z in v];g:=chosen[1]*xy[1]+chosen[2]*xy[2];fs:=[Evaluate(pol,g)]cat[xy[j]-Evaluate(PK!cs[j],g):j in[1,2]];J:=ideal<OF|fs>;pp:=[Place(z[1]):z in Factorization(J)|z[2]gt 0];if #pp eq 0 then pp:=[Place(z[1]):z in Factorization(ideal<OI|fs>)|z[2]gt 0];end if;assert #pp eq 1;P:=pp[1];assert Degree(P)eq e and &and[z eq 0 or Valuation(z,P)gt 0:z in fs];return P,d div e;
end function;
function SparseDiv(v)
 D:=DivisorGroup(F)!0;for z in v do assert z[1]in[1..#fb];D+:=z[2]*Divisor(fb[z[1]]);end for;return D;
end function;
function SparseClass(v)
 return &+[G|z[2]*fc[z[1]]:z in v];
end function;
out:=[* *];W:=eval Read(root cat "fourmap_class_witnesses.m");wi:=0;direct:=0;
for ix in[1..#packet[5]]do im:=packet[5][ix];P,m:=FromKey(im[1]);e:=Degree(P);assert im[3]eq e;ii:=e le Maximum(a[3])select Index(fb,P)else 0;target:=q(Z!im[2]);
 if ii gt 0 then assert target eq m*(fc[ii]-e*fc[ia]);direct+:=1;used:=0;
 else wi+:=1;assert wi le #W;v:=W[wi];assert #v eq 7 and v[2]eq 0 and m*q(Z!v[3])eq target;E:=SparseDiv(v[4]);R:=SparseDiv(v[5]);AA:=SparseDiv(v[6]);assert Degree(E)eq 10 and Degree(R)eq 10 and AA eq A;assert &and[z[2]ge 0:z in v[4]]and &and[z[2]ge 0:z in v[5]];b:=FF(v[7]);assert b ne 0;IE,JE:=Ideals(E);IR,JR:=Ideals(R);ID,JD:=Ideals(Divisor(P)-e*A);assert IE/(IR*(b*OF))eq ID and JE/(JR*(b*OI))eq JD;assert m*(SparseClass(v[4])-SparseClass(v[5]))eq target;used:=wi;
 end if;Append(~out,[*ix,e,m,ii,used,Eltseq(target)*]);if ix mod 50 eq 0 then print "C56_NATIVE_IDEAL_PROGRESS",ix,wi,direct,Cputime(tt);end if;
end for;assert wi eq #W;Write(root cat "independent_native_fourmap_ideals.m",Sprint(out):Overwrite:=true);print "D31B_C56_NATIVE_IDEAL_COMPLETE",#out,wi,direct,Cputime(tt);quit;
