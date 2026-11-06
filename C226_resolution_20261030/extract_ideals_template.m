SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(5*10^9);tt:=Cputime();SetUseMontes(FldFun,true);k:=Rationals();PP:=PolynomialRing(k);B<t>:=FunctionField(k);R<Y>:=PolynomialRing(B);model:=eval Read("work/deg10/C226_resolution_20261030/Q_model_fast.m");
function EQ(v) return k!v[1]/k!v[2];end function;
function ER(v) return B!(PP![EQ(a):a in v[1]])/B!(PP![EQ(a):a in v[2]]);end function;
F<z>:=FunctionField(R![ER(v):v in model[1]]);
function EF(v) return F![ER(a):a in v];end function;
function SQ(v) return [Numerator(v),Denominator(v)];end function;
OF:=MaximalOrderFinite(F:Al:="Montes");OI:=MaximalOrderInfinite(F);b:=EF(model[2]);c:=EF(model[3]);data:=eval Read(Sprintf("work/deg10/C226_resolution_20261030/effective_ideals_%o.m",CLASS));IE:=ideal<OF|[EF(v):v in data[1]]>;JE:=ideal<OI|[EF(v):v in data[2]]>;assert IsIntegral(IE) and IsIntegral(JE);pls:=[];mults:=[];
for fa in [*Factorization(IE),Factorization(JE)*] do for r in fa do Append(~pls,Place(r[1]));Append(~mults,r[2]);end for;end for;
assert &+[Degree(pls[i])*mults[i]:i in [1..#pls]] eq 10;print "Q226_IDEAL_SUPPORT",CLASS,[Degree(p):p in pls],mults,Cputime(tt);records:=[* *];
for i in [1..#pls] do
 p:=pls[i];K:=ResidueClassField(p);bk:=Evaluate(b,p);ck:=Evaluate(c,p);E:=EllipticCurve([K|0,ck,0,bk*(1-bk-ck),0]);P:=E![0,0];Q:=E![bk,bk];assert Discriminant(E) ne 0 and 2*P eq E!0 and P ne E!0 and 26*Q eq E!0 and 13*Q ne E!0 and 2*Q ne E!0 and P ne 13*Q;
 kp:=Type(K) eq FldRat select [k|0,1] else Coefficients(DefiningPolynomial(K));Append(~records,[*[SQ(v):v in kp],[SQ(v):v in Eltseq(K!bk)],[SQ(v):v in Eltseq(K!ck)],mults[i]*]);
end for;
Write(Sprintf("work/deg10/C226_resolution_20261030/exact_points_%o.m",CLASS),Sprint(records):Overwrite:=true);print "Q226_IDEAL_EXTRACTION_COMPLETE",CLASS,Cputime(tt);quit;
