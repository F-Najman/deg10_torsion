SetColumns(0);SetQuitOnError(true);SetMemoryLimit(2*10^9);tt:=Cputime();out:=Open("work/degree9/cyclic/continuation_20261026/C66_quotient17/source33/odd_order_bound33.jsonseq","w");bound:=0;expected:=5*42373650;
for q in [5,7,13,17,19,23,29,31,37,41,43]do
 f:=FrobeniusPolynomial(JOne(33),q);assert Degree(f)eq 42 and Coefficient(f,42)eq 1;ord:=Integers()!Evaluate(f,1);bound:=GCD(bound,ord);assert bound mod expected eq 0;
 fprintf out,"{\"prime\":%o,\"Frobenius_coefficients\":%o,\"order\":%o,\"gcd\":%o}\n",q,Coefficients(f),ord,bound;Flush(out);printf "NONIC33_ODD_ORDER_BOUND q=%o bound=%o ratio=%o time=%o\n",q,bound,bound div expected,Cputime(tt);if bound eq expected then break;end if;
end for;delete out;printf "NONIC33_ODD_ORDER_BOUND_COMPLETE bound=%o target=%o equality=%o time=%o\n",bound,expected,bound eq expected,Cputime(tt);quit;
