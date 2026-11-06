SetColumns(0);SetQuitOnError(true);SetNthreads(1);tt:=Cputime();AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
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

AF,af:=AlgorithmicFunctionField(FF);assert Degree(AF) eq 6;assert af(FF.1) eq AF!(BaseField(AF).1) and af(FF.2) eq AF.1;
Write("work/deg10/continuation_20261029/M6N12_char0_coordinate_identity.txt",Sprint([*Degree(AF),true,true*]));print "DECIC612_CHAR0_COORDINATE_IDENTITY_COMPLETE",Cputime(tt);quit;
