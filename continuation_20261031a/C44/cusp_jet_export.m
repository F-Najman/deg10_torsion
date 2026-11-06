SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(12*10^9);tt:=Cputime();
load "work/deg10/continuation_20261031a/C44/inputs/canonical_data.m";
k:=GF(3);K<t>:=FunctionField(k);PP<z>:=PolynomialRing(K);PR:=PolynomialRing(k);
function DR(a)return K!(PR!a[1])/(K!(PR!a[2]));end function;
f:=PP![DR(a):a in AFPolynomialData];F<u>:=FunctionField(f);assert Genus(F) eq 36;
function DF(a)return F![DR(v):v in a];end function;
xy:=[DF(a):a in ModularCoordinateData];x,y:=Explode(xy);
rr:=(x^2*y-x*y+y-1)/(x^2*y-x);ss:=(x*y-y+1)/(x*y);bb:=rr*ss*(rr-1);cc:=ss*(rr-1);
Eu:=EllipticCurve([1-cc,-bb,-bb,F!0,F!0]);Pu:=Eu![0,0];assert 44*Pu eq Eu!0 and 22*Pu ne Eu!0 and 4*Pu ne Eu!0;
jF:=jInvariant(Eu);cs:=Poles(jF);assert Degree(jF) eq 720 and {*Degree(Q):Q in cs*} eq {*1^^15,2^^5,5^^3,10*};
Z:=[DF(a):a in CanonicalCoordinateData];assert #Z eq 36;
printf "D31A_C44_CUSP_MODEL_READY signature=%o time=%o\n",[Degree(Q):Q in cs],Cputime(tt);
output:="work/deg10/continuation_20261031a/C44/cusp_jets.m";
Write(output,Sprintf("Ds:=%o; Widths:=%o; BasisJets:=[* *];",[Degree(Q):Q in cs],[-Valuation(jF,Q):Q in cs]));
for i in [1..#cs] do
 Q:=cs[i];dd:=Degree(Q);m:=10 div dd;ord:=Min([Valuation(v,Q):v in Z]);prec:=m+Abs(ord)+12;
 repeat
  LS,ex:=Completion(F,Q:Precision:=prec);zz:=LS.1;vv:=[zz^(-ord)*ex(v):v in Z];
  okay:=&and[(IsZero(v) or Valuation(v) ge 0) and AbsolutePrecision(v) ge m:v in vv];
  if not okay then prec*:=2;end if;
 until okay;
 rows:=[];for v in vv do row:=[];for j in [0..m-1] do row cat:=Eltseq(Coefficient(v,j),k);end for;Append(~rows,row);end for;
 A:=Matrix(k,rows);assert Nrows(A) eq 36 and Ncols(A) eq m*dd;
 Write(output,Sprintf("Append(~BasisJets,Matrix(GF(3),36,%o,%o));",m*dd,Eltseq(A)));
 printf "D31A_C44_CUSP_JETS index=%o/%o degree=%o multiplicity=%o rank=%o time=%o\n",i,#cs,dd,m,Rank(A),Cputime(tt);
end for;
printf "D31A_C44_CUSP_JETS_COMPLETE places=%o time=%o\n",#cs,Cputime(tt);quit;
