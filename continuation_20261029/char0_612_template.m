SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(8*10^9);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
data:=eval Read("work/deg10/continuation_20261029/M6N12_D5_all_source_lifts.m");
R<T>:=PolynomialRing(Rationals());F<z>:=NumberField(T^2-T+1);
 A<u,v>:=AffineSpace(F,2);
 rawmodel:=Split(Read("work/deg10/shared/mdmagma/models_X1_m_n/X1_6_12.txt"));
 for line in rawmodel do
  val:=Split(Split(line,"=")[2],";")[1];
  if line[1] eq "X" then equation:=eval(val);end if;
  if line[1] eq "r" then r:=eval(val);end if;
  if line[1] eq "s" then s:=eval(val);end if;
  if line[1] eq "q" then q:=eval(val);end if;
  if line[1] eq "E" then es:=eval(val);end if;
 end for;
 C:=ProjectiveClosure(Curve(A,equation));FF:=FunctionField(C);
 X:=New(MDCrvMod11);X`M:=6;X`N:=12;
 X`_P:=[FF!q];X`_Q:=[FF!0,FF!0];X`_coordinates:=[FF.1,FF.2];
 _InitMDCrvMod(X,12,C,[FF!a:a in es]);

C:=Curve(X);FF:=FunctionField(C);print "M6N12_CHAR0_MODEL",Genus(X),Cputime(tt);
 E:=EllipticCurve([FF!a:a in X`_E]);Q:=E![0,0];assert 12*Q eq E!0;
 px:=FF!X`_P[1];ai:=aInvariants(E);R<Y>:=PolynomialRing(FF);
 ys:=Roots(Y^2+(ai[1]*px+ai[3])*Y-px^3-ai[2]*px^2-ai[4]*px-ai[5]);
 assert #ys eq 2;P:=E![px,ys[1][1]];

 if WeilPairing(P,2*Q,6) ne z then P:=-P;end if;
 assert WeilPairing(P,2*Q,6) eq z;
 print "M6N12_CHAR0_WEIL",Cputime(tt);
 fs:=[(P+5*Q)[1]-P[1],(2*P+5*Q)[1]-P[1],(P+2*Q)[1]];
 assert &and[f ne 0:f in fs];print "M6N12_CHAR0_THREE_UNITS",Cputime(tt);
 cs:=Cusps(X);groups:=AssociativeArray();
 for c in cs do
  s:=[Valuation(f,c):f in fs];
  if IsDefined(groups,s) then Append(~groups[s],c);else groups[s]:=[c];end if;
 end for;
 cols:=[34,58,27];profiles:=[[sig[j]:j in cols]:sig in data[1]];
 assert #Seqset(profiles) eq 24 and Seqset(profiles) eq Keys(groups);
 blocks:=[&+[Divisor(c):c in groups[sig]]:sig in profiles];ds:=[Degree(D):D in blocks];
 print "M6N12_CHAR0_CUSP_BLOCKS",#blocks,ds,Cputime(tt);
 sigs:=data[1];
 c0:=Support(blocks[Index(ds,1)])[1];

results:=[* *];
for index in CHOSEN do
 coeff:=data[3][index][2];L:=10*Divisor(c0)+&+[coeff[i]*(blocks[i]-ds[i]*Divisor(c0)):i in [1..24]];assert Degree(L) eq 10;
 V,mp:=RiemannRochSpace(L);assert Dimension(V) eq 2;ff:=mp(V.1)/mp(V.2);Write(Sprintf("work/deg10/continuation_20261029/char0_612_class%o_function.m",index),Sprint(ff));
 ps,es:=Support(Divisor(ff));poles:=&+[-es[k]*Divisor(ps[k]):k in [1..#ps]|es[k] lt 0];degree:=Degree(poles);assert degree le 10;
 print "DECIC612_CHAR0_PENCIL",index,degree,Cputime(tt);ramdeg:=[];rammult:=[];quintic:=[* *];
 if degree eq 10 then
  ram:=Divisor(Differential(ff))+2*poles;assert IsEffective(ram) and Degree(ram) eq 36;rp,rm:=Support(ram);ramdeg:=[Degree(P):P in rp];rammult:=rm;
  print "DECIC612_CHAR0_RAMIFICATION",index,ramdeg,rammult,Cputime(tt);
  for k in [1..#rp] do
   P:=rp[k];if Degree(P) ne 5 then continue;end if;
   if Valuation(ff,P) lt 0 then isF:=true;value:=F!0;fiber:=poles;
   else val:=Evaluate(ff,P);isF,value:=IsCoercible(F,val);if isF then fiber:=poles+Divisor(ff-value);end if;end if;
   if isF then
    assert fiber eq 2*Divisor(P);Ej:=EllipticCurve(X,P);j:=jInvariant(Ej);print "DECIC612_CHAR0_QUINTIC_POINT_FOUND",index,P,j,Cputime(tt);
    Append(~quintic,[*Sprint(P),Sprint(aInvariants(Ej)),Sprint(j)*]);
    Write(Sprintf("work/deg10/continuation_20261029/char0_612_class%o_quintic_point.txt",index),Sprint([*P,aInvariants(Ej),j*]));
   end if;
  end for;
 end if;
 Append(~results,[*index,degree,ramdeg,rammult,quintic*]);
end for;
Write("OUTFILE",Sprint(results));print "DECIC612_CHAR0_PENCILS_COMPLETE",results,Cputime(tt);quit;
