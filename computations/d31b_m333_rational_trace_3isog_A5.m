SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(6*10^9);Alarm(850);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");Attach("work/deg10/continuation_20261031b/M3N33/inputs/fast_hecke.m");
n:=33;root:="work/deg10/continuation_20261031b/M3N33/";k:=GF(2);PP:=PolynomialRing(GF(2));X:=MDX1(n,k);F0,mp:=AlgorithmicFunctionField(FunctionField(Curve(X)));F:=RationalExtensionRepresentation(F0);coords:=[F!mp(v):v in X`_coordinates_xy];OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);print "RELATIVE5_TARGET_MODEL",n,Genus(F),Degree(F),Cputime(tt);
function Enc(a) return Integers()!(k!a);end function;
function ER(v)return [[Enc(c):c in Coefficients(Numerator(v))],[Enc(c):c in Coefficients(Denominator(v))]];end function;
function EF(v)return [ER(c):c in Eltseq(F!v)];end function;
function PD(p)
 O:=Order(Ideal(p));aa,bb:=TwoGenerators(p);ts:=[aa,bb];assert ideal<O|ts> eq Ideal(p);return [*Degree(p),O eq OF select 1 else 0,[EF(v):v in ts]*];end function;
function Get(E,P)
 values:=TateNormalForm_xy(E,P);K:=Parent(values[1]);Embed(k,K);e:=LCM([Degree(MinimalPolynomial(v,k)):v in values]);found:=false;
 for bits in [[k!1,k!0],[k!0,k!1]] cat [[k!1,a]:a in k|a ne 0] do
  gen:=&+[K!bits[i]*values[i]:i in[1,2]];pol:=MinimalPolynomial(gen,k);if Degree(pol) eq e then chosen:=bits;found:=true;break;end if;end for;assert found;
 M:=Matrix(k,[Eltseq(gen^i,k):i in[0..e-1]]);assert Rank(M) eq e;cs:=[Eltseq(Solution(M,Vector(k,Eltseq(v,k)))):v in values];assert &and[&+[K!cs[j][i+1]*gen^i:i in[0..e-1]] eq values[j]:j in[1,2]];
 gg:=&+[F!chosen[i]*coords[i]:i in[1,2]];fs:=[Evaluate(pol,gg)] cat [coords[j]-&+[F!cs[j][i+1]*gg^i:i in[0..e-1]]:j in[1,2]];I:=ideal<OF|fs>;fa:=Factorization(I);ps:=[Place(z[1]):z in fa|z[2] gt 0];if #ps eq 0 then fa:=Factorization(ideal<OI|fs>);ps:=[Place(z[1]):z in fa|z[2] gt 0];end if;assert #ps eq 1 and Degree(ps[1]) eq e;assert &and[f eq 0 or Valuation(f,ps[1]) gt 0:f in fs];return ps[1],[[Integers()!v:v in Eltseq(t,GF(2))]:t in values];
end function;
data:=eval Read(root cat "physical_source5.m");K<z>:=ext<GF(2)|PP!data[1]>;tests:=[* *];
Write(root cat "rational_trace_model.m",Sprint([*Coefficients(DefiningPolynomial(k)),[ER(v):v in Coefficients(DefiningPolynomial(F))],[EF(v):v in coords]*]):Overwrite:=true);
for i in[1,6,10] do
 rec:=data[6][i];ky:=rec[4];b:=Evaluate(PP!ky[1..10],z);c:=Evaluate(PP!ky[11..20],z);E:=EllipticCurve([K|1-c,-b,-b,0,0]);Q:=E![0,0];assert Order(Q) eq n;Eorig:=E;Qorig:=Q;RX:=PolynomialRing(K);r3:=Roots(DivisionPolynomial(Eorig,3));assert #r3 eq 4;
 for rr3 in r3 do
 kerx:=rr3[1];E,phi3:=IsogenyFromKernel(Eorig,RX.1-kerx);Q:=phi3(Qorig);if Order(Q) ne n then assert Order(Q) eq 11;continue;end if;
 ps,xy:=Get(E,Q);pd,xyd:=Get(E,5*Q);ns:=10 div Degree(ps);nd:=50 div Degree(pd);assert ns*Degree(ps) eq 10 and nd*Degree(pd) eq 50;D:=-ns*Divisor(ps)-nd*Divisor(pd);hec:=[* *];fast:=MDIsogeniesFast(E,Q,5);assert &+[v[3]:v in fast] eq 6;
 for row in fast do ET,QT,e:=Explode(row);pt,xyt:=Get(ET,QT);nn:=10*e div Degree(pt);assert nn*Degree(pt) eq 10*e;D+:=nn*Divisor(pt);Append(~hec,[*e,nn,PD(pt),xyt,[Integers()!v:v in Coefficients(DefiningPolynomial(BaseRing(ET),GF(2)))]*]);end for;
 assert Degree(D) eq 0;h0:=#Basis(D:Reduction:=false,Simplification:="None");Append(~tests,[*i,Eltseq(kerx,GF(2)),[[Integers()!a:a in Eltseq(v,GF(2))]:v in aInvariants(E)],[[Integers()!a:a in Eltseq(v,GF(2))]:v in Eltseq(Q)],ky,PD(ps),ns,PD(pd),nd,xy,xyd,hec,h0,rec[7]*]);Write(root cat "rational_trace_3isog_A5_progress.m",Sprint(tests):Overwrite:=true);print "RELATIVE5_RATIONAL_TRACE_3ISOG_A5_RR",n,i,Eltseq(kerx,GF(2)),h0,Cputime(tt);
 end for;
end for;
Write(root cat "rational_trace_3isog_A5_RR.m",Sprint(tests):Overwrite:=true);print "D31B_RELATIVE5_RATIONAL_TRACE_A5_COMPLETE",n,#tests,[v[13]:v in tests],Cputime(tt);quit;
