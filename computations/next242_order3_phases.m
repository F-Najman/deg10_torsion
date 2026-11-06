SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(10^9);k:=GF(5);R:=PolynomialRing(k);K<a>:=ext<k|R![3,1,0,2,0,0,0,0,0,1]>;
rows:=[*[1777,1,0,0,0,0,2,0,2,1,4,4,4,3,1,0,0,4,2,2],[1780,2,0,0,1,0,1,3,4,1,2,3,4,1,0,4,1,1,4,4],[1783,2,0,0,0,1,0,1,0,2,3,4,2,3,2,3,3,0,2,1]*];out:=[];
function BC(E,P,Q) x:=Q[1]-P[1];y:=Q[2];return [x^3/y^2,(aInvariants(E)[2]+3*P[1])*x^2/y^2];end function;
function Act(P,Q,n) if n eq 1 then return 21*Q,P+25*Q;else return P+21*Q,P+4*Q;end if;end function;
for row in rows do ci,n:=Explode(row[1..2]);b:=Evaluate(R!row[3..11],a);c:=Evaluate(R!row[12..20],a);E:=EllipticCurve([K|0,c,0,b*(1-b-c),0]);P:=E![0,0];Q:=E![b,b];assert 42*Q eq E!0 and 21*Q ne E!0 and 14*Q ne E!0 and 6*Q ne E!0 and P ne 21*Q;P1,Q1:=Act(P,Q,n);u,v:=Act(P1,Q1,n);u,v:=Act(u,v,n);assert u eq P and v eq Q;bc:=BC(E,P1,Q1);phases:=[j:j in [0..8]|bc eq [b^(5^j),c^(5^j)]];assert #phases eq 1 and phases[1] in [3,6];Append(~out,[ci,n,phases[1]]);end for;
Write("work/deg10/next_cases_20261030/C242/degree9_order3_native_phases.m",Sprint(out):Overwrite:=true);print "C242_THREE_ORDER3_NATIVE_PHASES_PASS",out;quit;
