SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(3*10^9);t0:=Cputime();
J:=JH(54,[19]);MM:=ModularSymbols(J);assert #MM eq 1;M:=MM[1];assert Dimension(M) eq 20;
B:=Matrix(Rationals(),BasisMatrix(IntegralHomology(J)));assert Nrows(B) eq 20 and Ncols(B) eq 20;
S:=B*Matrix(Rationals(),StarInvolution(M))*B^-1;assert &and[Denominator(a) eq 1:a in Eltseq(S)];S:=Matrix(Integers(),S);assert S^2 eq 1;
I:=IdentityMatrix(Integers(),20);ops:=[* *];rowsplus:=Transpose(S-I);rowsminus:=Transpose(S+I);pack:=[* *];
for q in [5,7,11,13,17,19,23,31] do
 T:=Matrix(Rationals(),HeckeOperator(M,q));DD:=(T^2-Matrix(Rationals(),HeckeOperator(M,q^2)))/q;
 T:=B*T*B^-1;DD:=B*DD*B^-1;
 assert &and[Denominator(a) eq 1:a in Eltseq(T) cat Eltseq(DD)];T:=Matrix(Integers(),T);DD:=Matrix(Integers(),DD);
 assert S*T eq T*S and S*DD eq DD*S;
 ep:=T-I-q*DD;em:=T+(q mod 3 eq 1 select -1 else 1)*(I+q*DD);
 rowsplus:=VerticalJoin(rowsplus,Transpose(ep));rowsminus:=VerticalJoin(rowsminus,Transpose(em));
 Sp:=SmithForm(rowsplus);Sm:=SmithForm(rowsminus);dp:=[Sp[i,i]:i in [1..20]];dm:=[Sm[i,i]:i in [1..20]];assert &and[x ne 0:x in dp cat dm];np:=&*dp;nm:=&*dm;
 Append(~pack,[*q,Eltseq(T),Eltseq(DD),dp,dm*]);
 print "LATE318_HOMOLOGY_TORSION_BOUND",q,np,nm,Valuation(np,3)+Valuation(nm,3),Cputime(t0);
 Write("work/deg10/continuation_20261029b/homology318.m",Sprint([*Eltseq(B),Eltseq(S),pack*]):Overwrite:=true);
 if Valuation(np,3)+Valuation(nm,3) eq 10 then break;end if;
end for;
print "LATE318_HOMOLOGY_COMPLETE",Cputime(t0);quit;
