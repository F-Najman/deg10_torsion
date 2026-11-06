SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(8*10^9);Alarm(240);tt:=Cputime();root:="work/deg10/continuation_20261031b/C56/";stems:=["cover_v4_A5_s0","cover_v4_A5_s1","cover_v4_A5_s2","cover_v4_A11_s0"];md:=eval Read(root cat stems[1] cat "_model.m");k:=GF(3);PK:=PolynomialRing(k);K<t>:=FunctionField(k);R<X>:=PolynomialRing(K);
function DR(v)return K!(PK!v[1])/(K!(PK!v[2]));end function;
F<u>:=FunctionField(R![DR(v):v in md[1]]);assert Genus(F) eq 25;OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
function DF(v)return F![DR(a):a in v];end function;
function BaseP(v)OO:=v[3] eq 1 select OF else OI;gs:=[DF(a):a in v[4]];II:=ideal<OO|gs>;assert #gs eq 2 and IsPrime(II);p:=Place(II);assert Degree(p) eq v[2];return p;end function;
function NF(v)assert v[1] mod v[2] eq 0;return v[1] div v[2];end function;
function ER(v)return [[Integers()!a:a in Coefficients(Numerator(v))],[Integers()!a:a in Coefficients(Denominator(v))]];end function;
B:=DF(md[2][1]);branch:=[p:p in Zeros(B) cat Poles(B)|IsOdd(Valuation(B,p))];assert &+[Degree(p):p in branch] eq 24;out:=[* *];
for stem in stems do
 dat:=eval Read(root cat stem cat "_operator_progress.m");rows,cas,inds:=Explode(dat);D:=DivisorGroup(F)!0;for i in inds do row:=rows[i];Q:=row[7];D-:=NF(row[4])*Divisor(BaseP(row[4]));D-:=Q*NF(row[6])*Divisor(BaseP(row[6]));for h in row[8] do D+:=h[2]*NF(h[3])*Divisor(BaseP(h[3]));end for;end for;assert Degree(D) eq 0;rr:=Basis(D);h0:=#rr;print "C56_BASE_NORM_RR",cas[1],h0,Cputime(tt);signs:=[];valid:=[];gg:=[];
 if h0 eq 1 then g:=rr[1];assert Divisor(g)+D eq 0;for p in branch do assert Valuation(g,p) eq 0;v:=Evaluate(g,p);chi:=v^((3^Degree(p)-1) div 2);assert chi eq 1 or chi eq -1;Append(~signs,[Degree(p),chi eq 1 select 1 else -1]);end for;valid:=[c:c in[1,2]|&and[s[2]*(c eq 1 or IsEven(s[1]) select 1 else -1) eq 1:s in signs]];gg:=[ER(v):v in Eltseq(g)];end if;
 Append(~out,[*cas[1],Q,h0,valid,signs,gg*]);print "C56_BRANCH_NORM_OBSTRUCTION",cas[1],valid,signs,Cputime(tt);Write(root cat "branch_norm_probe.m",Sprint(out):Overwrite:=true);
end for;print "D31B_C56_BRANCH_NORM_PROBE_COMPLETE",Cputime(tt);quit;
