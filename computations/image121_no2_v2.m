/* Isogeny-invariant reduction orders on A=(diamond16-1)J1(115).
   Oddness of #A(Fq) excludes rational2-primary torsion independently
   of any choice of integral homology lattice. */
SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(4*10^9);tt:=Cputime();
N:=121;aa:=2;D:=FullDirichletGroup(N);chiA:=[];genusA:=0;
for chi in GaloisConjugacyRepresentatives(D) do
 if IsOdd(chi) or Evaluate(chi,aa) eq 1 then continue;end if;
 chimin:=MinimalBaseRingCharacter(chi);Append(~chiA,chimin);
 S:=CuspidalSubspace(ModularSymbols(chimin,2,1));genusA+:=Dimension(S)*EulerPhi(Order(chi));
end for;
printf "DECIC_121_DIAMOND_IMAGE dimension=%o character_orders=%o time=%o\n",genusA,[Order(c):c in chiA],Cputime(tt);
for q in [3] do
 allodd:=true;
 for chi in chiA do
  S:=CuspidalSubspace(ModularSymbols(chi,2,1));d:=Dimension(S);Tq:=HeckeOperator(S,q);K:=BaseRing(Tq);
  A:=(1+q*(K!Evaluate(chi,q)))*IdentityMatrix(K,d)-Tq;
  detA:=Determinant(A);z:=Norm(detA);assert Denominator(z) eq 1 and z ne 0;
  odd:=IsOdd(Integers()!z);allodd:=allodd and odd;
  rr:=func<el|[[Numerator(c),Denominator(c)]:c in Eltseq(K!el)]>;
  pol:=[[Numerator(c),Denominator(c)]:c in Coefficients(DefiningPolynomial(K))];
  file:=Sprintf("work/deg10/exclusions/image121_chi%o.json",Order(chi));
  Write(file,Sprintf("{\"N\":121,\"q\":3,\"chi_order\":%o,\"d\":%o,\"field_polynomial\":%o,\"A\":%o,\"norm\":%o,\"norm_power\":%o}",Order(chi),d,pol,[rr(el):el in Eltseq(A)],z,Degree(K) div EulerPhi(Order(chi))));
  printf "DECIC_121_IMAGE_REDUCTION q=%o character_order=%o relative_dimension=%o coefficient_degree=%o norm_power=%o norm=%o odd=%o time=%o\n",q,Order(chi),d,Degree(K),Degree(K) div EulerPhi(Order(chi)),z,odd,Cputime(tt);
 end for;
 assert allodd;
 if allodd then
  printf "DECIC_121_IMAGE_NO_TWO_PRIMARY_PASS q=%o complete=true time=%o\n",q,Cputime(tt);break;
 end if;
end for;
print "DECIC_121_IMAGE_EVEN_ORDERS_CERTIFIED_COMPLETED";quit;
