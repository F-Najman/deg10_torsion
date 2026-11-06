SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(2*10^9);tt:=Cputime();
load "work/deg10/continuation_20261030c/C91_CM_control/inputs/records91_8.m";k:=GF(2);PK:=PolynomialRing(k);K<z>:=ext<k|PK!FieldPolynomial>;RX<X>:=PolynomialRing(K);zetas:=[v[1]:v in Roots(X^2+X+1)];assert #zetas eq 2;out:=[* *];
for i in[1..#Records] do rec:=Records[i];E:=EllipticCurve([Evaluate(PK!v,z):v in rec[2]]);P:=E![Evaluate(PK!v,z):v in rec[3]];aa:=aInvariants(E);assert aa[1] eq 0 and aa[2] eq 0 and aa[4] eq 0 and aa[5] eq 0 and Order(P) eq 91;
 zz:=[u:u in zetas|E![u*P[1],P[2]] eq 16*P];assert #zz eq 1;u:=zz[1];P7:=13*P;S13:=7*P;assert Order(P7) eq 7 and Order(S13) eq 13;assert E![u*P7[1],P7[2]] eq 2*P7;assert E![u*S13[1],S13[2]] eq 3*S13;
 ker:=&*[X-(j*S13)[1]:j in[1..6]];autker:=&*[X-u*(j*S13)[1]:j in[1..6]];assert autker eq ker;Append(~out,[*i,[Integers()!a:a in Eltseq(u,k)],16,2,3,true*]);end for;
Write("work/deg10/continuation_20261030c/C91_CM_control/CM_level7_binding.m",Sprint(out):Overwrite:=true);print "C91_CM_AUTOMORPHISM_POINT7_KERNEL13_BINDING_COMPLETE",#out,Cputime(tt);quit;
