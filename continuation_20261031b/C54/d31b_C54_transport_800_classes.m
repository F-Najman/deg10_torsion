SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(6*10^9);Alarm(1000);tt:=Cputime();AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");root:="work/deg10/continuation_20261031b/C54/";a:=eval Read(root cat "dyadic_rigidity_v2.m");X:=MDX1(27,GF(5));C:=Curve(X);F0,mp:=AlgorithmicFunctionField(FunctionField(C));F:=RationalExtensionRepresentation(F0);assert Genus(F)eq 13;cs:=Cusps(X);sigs:=[CuspSignature(X,p):p in cs];cp:=[Support(FunctionFieldDivisor(Divisor(cs[Index(sigs,sg)])))[1]:sg in a[3]];assert [Degree(p):p in cp]eq a[4];base:=Index(a[4],1);A:=Divisor(cp[base]);G:=AbelianGroup(a[1]);Z:=FreeAbelianGroup(#cp);hh:=hom<Z->G|[G!v:v in a[6]]>;assert #Image(hh)eq 3*3*52497;
classes:=Sort(Setseq({r[7]:r in a[8]}));assert #classes eq 800;vv:=[Eltseq((G!v)@@hh):v in classes];mx:=[Max([Abs(v[j]):v in vv]):j in[1..#cp]];powers:=AssociativeArray();
function RD(D)
 E,r,aa,h:=Reduction(D,A);assert aa eq A;B:=E+r*A;assert Degree(B)eq 0;return B;
end function;
for j in[1..#cp]do powers[j]:=[];D:=Divisor(cp[j])-Degree(cp[j])*A;n:=1;while n le mx[j]do Append(~powers[j],D);D:=RD(2*D);n*:=2;end while;end for;
print "D31B_C54_800_TRANSPORT_FRAME_READY",#classes,mx,Cputime(tt);OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
function ER(v)return[[Integers()!z:z in Coefficients(Numerator(v))],[Integers()!z:z in Coefficients(Denominator(v))]];end function;
function EF(v)return[ER(z):z in Eltseq(F!v)];end function;
function PD(p)return[*Degree(p),Order(Ideal(p))eq OF select 1 else 0,[EF(z):z in Basis(Ideal(p))]*];end function;
function DD(D)return[[*PD(z[1]),z[2]*]:z in Decomposition(D)];end function;
Write(root cat "transport_800_model.m",Sprint([*[ER(v):v in Coefficients(DefiningPolynomial(F))],[[EF(F!mp(z)):z in X`_coordinates_xy],[EF(F!mp(z)):z in X`_coordinates]],DD(A),[DD(Divisor(p)-Degree(p)*A):p in cp]*]):Overwrite:=true);
obs:=[* *];left:=[* *];for ci in[1..#classes]do D:=DivisorGroup(F)!0;v:=vv[ci];for j in[1..#cp]do n:=Abs(v[j]);s:=Sign(v[j]);ix:=1;while n gt 0 do if n mod 2 eq 1 then D:=RD(D+s*powers[j][ix]);end if;n div:=2;ix+:=1;end while;end for;D+:=10*A;bs:=Basis(D);h:=#bs;row:=[*ci,classes[ci],v,DD(D),h,[EF(b):b in bs]*];Append(~obs,row);if h gt 0 then Append(~left,ci);end if;if ci mod 25 eq 0 then print "D31B_C54_800_RR_PROGRESS",ci,#left,Cputime(tt);Write(root cat "transport_800_progress.m",Sprint([*obs,left*]):Overwrite:=true);end if;end for;
Write(root cat "transport_800_classes.m",Sprint([*obs,left*]):Overwrite:=true);print "D31B_C54_800_RR_COMPLETE",#obs,#left,Cputime(tt);quit;
