SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(10^9);root:="work/deg10/next_cases_20261030/C242/";aa:=eval Read(root cat "direct_hecke17/input_cuspidal.m");r:=[r:r in aa[2]|r[1] eq 26][1][2][1];d,ky:=Explode(r);assert d eq 6;pp:=[r[2]:r in aa[1]|r[1] eq d][1];k:=GF(5);R:=PolynomialRing(k);K<a>:=ext<k|R!pp>;b:=Evaluate(R!ky[1..d],a);c:=Evaluate(R!ky[d+1..2*d],a);E:=EllipticCurve([K|0,c,0,b*(1-b-c),0]);P:=E![0,0];Q:=E![b,b];
function BC(E,P,Q) x:=Q[1]-P[1];y:=Q[2];return [x^3/y^2,(aInvariants(E)[2]+3*P[1])*x^2/y^2];end function;
rows:=[* *];for no in [1,2] do if no eq 1 then PP:=21*Q;QQ:=P+25*Q;else PP:=P+21*Q;QQ:=P+4*Q;end if;bc:=BC(E,PP,QQ);ph:=[i:i in [0..5]|bc eq [b^(5^i),c^(5^i)]];assert #ph eq 1;Append(~rows,[no,ph[1]]);end for;
// The reflection interchanges the two diagonal order3 subgroups.
PP:=21*Q;QQ:=P+22*Q;bc:=BC(E,PP,QQ);Write(root cat "case26_diagonal_phases.m",Sprint([*rows,bc eq [b,c]*]):Overwrite:=true);print "C242_CASE26_DIAGONAL_PHASES",rows;quit;
