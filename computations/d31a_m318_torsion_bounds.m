SetColumns(0);SetNthreads(1);SetQuitOnError(true);SetMemoryLimit(4*10^9);Alarm(200);tt:=Cputime();J:=JH(54,[19]);assert Dimension(J) eq 10;rows:=[];bound:=0;
for p in[5,7,11,13,19,31] do
 f:=FrobeniusPolynomial(J,p);assert Degree(f) eq 20;e:=p mod 3 eq 1 select 1 else 2;n:=Integers()!(Evaluate(f,1)*(e eq 1 select 1 else Evaluate(f,-1)));bound:=GCD(bound,n);Append(~rows,[p,e,n,Valuation(n,3),bound]);print "M318_TORSION_ORDER_BOUND",p,e,n,Valuation(n,3),bound,Cputime(tt);
end for;
Write("work/deg10/continuation_20261031a/M3N18/torsion_order_bounds.m",Sprint(rows):Overwrite:=true);print "M318_TORSION_ORDER_BOUNDS_COMPLETE",bound,Cputime(tt);quit;
