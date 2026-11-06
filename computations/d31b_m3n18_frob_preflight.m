SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(4*10^9);Alarm(450);tt:=Cputime();
J:=JH(54,[19]);g:=Dimension(J);print "RELATIVE_PREFLIGHT_GENUS",18,g,Cputime(tt);rows:=[* *];
for p in [2,7,13] do
 if 18 mod p eq 0 then continue;end if;
 f:=FrobeniusPolynomial(J,p);assert Degree(f) eq 2*g;tr:=[];e:=p mod 3 eq 1 select 1 else 2;q:=p^e;
 for i in[1..5*e] do
 a:=-i*Coefficient(f,2*g-i);for j in[1..i-1] do a-:=Coefficient(f,2*g-j)*tr[i-j];end for;Append(~tr,Integers()!a);end for;
 cl:=[];for d in[1..5] do a:=q^d+1-tr[e*d];for j in Divisors(d) do if j lt d then a-:=j*cl[j];end if;end for;assert a mod d eq 0;Append(~cl,a div d);end for;Append(~rows,[*p,e,cl,Coefficients(f)*]);print "RELATIVE_PREFLIGHT_CLOSED",p,e,cl,Cputime(tt);end for;
Write("work/deg10/continuation_20261031b/M3N18/frob_preflight.m",Sprint([*g,rows*]):Overwrite:=true);print "D31B_M3N18_FROB_PREFLIGHT_COMPLETE",Cputime(tt);quit;
