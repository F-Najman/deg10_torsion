SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(6*10^9);Alarm(400);tt:=Cputime();root:="work/deg10/continuation_20261031b/C54/";k:=GF(5);Rx:=PolynomialRing(k);R<x>:=FunctionField(k);S<y>:=PolynomialRing(R);a:=eval Read(root cat "transport_800_model.m");rows:=eval Read(root cat "class_members_profile_input.m");
function DR(v)return Evaluate(Rx!v[1],x)/Evaluate(Rx!v[2],x);end function;
F<t>:=FunctionField(S![DR(v):v in a[1]]);assert Genus(F)eq 13;OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
function DF(v)return F![DR(z):z in v];end function;
function DP(v)O:=v[2]eq 1 select OF else OI;I:=ideal<O|[DF(z):z in v[3]]>;p:=Place(I);assert Degree(p)eq v[1];return p;end function;
function DD(v)return &+[z[2]*Divisor(DP(z[1])):z in v];end function;
b,c:=Explode([DF(v):v in a[2][2]]);j:=jInvariant(EllipticCurve([F|1-c,-b,-b,0,0]));
function Key(vals)
 K:=Parent(vals[1]);d:=LCM([Degree(MinimalPolynomial(v,k)):v in vals]);chosen:=[];for z in[[1,0],[0,1]]cat[[1,n]:n in[1..4]]do gen:=z[1]*vals[1]+z[2]*vals[2];pol:=MinimalPolynomial(gen,k);if Degree(pol)eq d then chosen:=z;break;end if;end for;assert #chosen eq 2;M:=Matrix(k,[Eltseq(gen^i,k):i in[0..d-1]]);cs:=[Eltseq(Solution(M,Vector(k,Eltseq(v,k)))):v in vals];return[d]cat chosen cat[Integers()!v:v in Coefficients(pol)]cat &cat[[Integers()!v:v in r]:r in cs];end function;
points:=[* *];pointidx:=AssociativeArray();members:=[* *];pure:=0;outid:=0;
for row in rows do ci:=row[1];D:=DD(row[4]);assert Degree(D)eq 10;bs:=[DF(v):v in row[6]];assert #bs eq row[5];h:=#bs;
 for pivot in[1..h]do for vv in CartesianPower(k,h-pivot)do v:=[k|0:i in[1..pivot-1]]cat[k|1]cat[vv[i]:i in[1..h-pivot]];fn:=&+[v[i]*bs[i]:i in[1..h]];E:=D+Divisor(fn);assert Degree(E)eq 10 and &and[Valuation(E,p)ge 0:p in Support(E)];nc:=[];cd:=0;
  for p in Support(E)do mu:=Valuation(E,p);if Valuation(j,p)lt 0 then cd+:=mu*Degree(p);else vals:=[Evaluate(z,p):z in[b,c]];ky:=Key(vals);assert ky[1]eq Degree(p);if not IsDefined(pointidx,ky)then pointidx[ky]:=#points+1;Append(~points,[*#points+1,Degree(p),ky*]);end if;Append(~nc,[pointidx[ky],mu]);end if;end for;
  if #nc eq 0 then pure+:=1;else outid+:=1;Append(~members,[*outid,ci,row[2],nc,cd*]);end if;
 end for;end for;
 if ci mod 10 eq 0 then print "D31B_C54_RR_MEMBERS_PROGRESS",ci,#points,#members,pure,Cputime(tt);end if;
end for;
Write(root cat "class_members_profile.m",Sprint([*points,members,pure*]):Overwrite:=true);print "D31B_C54_RR_MEMBERS_PROFILE_COMPLETE",#rows,#points,#members,pure,Cputime(tt);quit;
