SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(6*10^9);
N:=1147;
// Prefix from bounded job log SHA256 1f1e43e69bc9c841cde0b5098ceb3961145e43ce71e08c258fb42348036236e3; only its completed character components are reused.
prefix:=[[1, 1, 99, 56, 1], [15, 31, 100, 100, 8], [5, 31, 98, 98, 4], [3, 31, 100, 100, 2], [180, 1147, 99, 99, 48], [180, 1147, 100, 100, 48], [36, 1147, 99, 99, 12], [180, 1147, 99, 99, 48], [36, 1147, 100, 100, 12], [36, 1147, 99, 99, 12], [18, 37, 100, 100, 6], [90, 1147, 99, 99, 24], [90, 1147, 99, 99, 24], [90, 1147, 100, 100, 24], [18, 1147, 99, 99, 6], [18, 1147, 99, 99, 6], [60, 1147, 100, 100, 16], [60, 1147, 98, 98, 16], [12, 1147, 100, 100, 4], [60, 1147, 100, 100, 16], [12, 1147, 98, 98, 4], [12, 1147, 100, 100, 4], [9, 37, 100, 100, 6], [45, 1147, 99, 99, 24], [45, 1147, 99, 99, 24], [45, 1147, 100, 100, 24], [9, 1147, 99, 99, 6], [9, 1147, 99, 99, 6], [6, 37, 98, 98, 2], [30, 1147, 100, 100, 8], [30, 1147, 100, 100, 8]];
D:=FullDirichletGroup(N);bad:=[];allchi:=[];total:=0;chiindex:=0;
for chi in GaloisConjugacyRepresentatives(D) do
 if IsOdd(chi) then continue;end if;chiindex+:=1;
 if chiindex le #prefix then
  r:=prefix[chiindex];assert Order(chi) eq r[1] and Conductor(chi) eq r[2];
  total+:=r[3]*EulerPhi(r[1]);Append(~allchi,chi);if r[3] ne r[4] then Append(~bad,chi);end if;
  printf "WINDING N=%o order=%o conductor=%o relative_dimension=%o winding_dimension=%o coefficient_degree=%o reused_completed_prefix=true\n",N,r[1],r[2],r[3],r[4],r[5];continue;
 end if;
 tt:=Cputime();
 chimin:=MinimalBaseRingCharacter(chi);
 S:=CuspidalSubspace(ModularSymbols(chimin,2,1));
 dd:=Dimension(S);
 printf "WINDING_SPACE_READY N=%o order=%o dimension=%o time=%o\n",N,Order(chi),dd,Cputime(tt);
 if Order(chi) ge 20 then
  vambient:=Vector(Eltseq(WindingElement(S)))*ProjectionMatrix(S);
  v:=Vector(Coordinates(VectorSpace(S),vambient));
  A:=HeckeOperator(S,2);
  printf "WINDING_T2_READY N=%o order=%o time=%o\n",N,Order(chi),Cputime(tt);
  K:=BaseRing(A);
  coeffs:=&cat[Eltseq(K!x):x in Eltseq(v) cat Eltseq(A)];
  pp:=0;
  for mult in [1..2000] do
   cand:=Order(chi)*mult+1;
   if IsPrime(cand) and &and[Denominator(c) mod cand ne 0:c in coeffs] then pp:=cand;break;end if;
  end for;
  assert pp gt 0;F:=GF(pp);
  printf "WINDING_GOOD_PRIME p=%o all_denominators_invertible=true\n",pp;
  roots:=Roots(PolynomialRing(F)!DefiningPolynomial(K));assert #roots gt 0;r:=roots[1][1];
  red:=func<x|&+[F!c*r^(i-1):i->c in Eltseq(K!x)]>;
  vm:=Vector(F,[red(x):x in Eltseq(v)]);
  Am:=Matrix(F,dd,dd,[red(x):x in Eltseq(A)]);
  vs:=[];u:=vm;
  for k in [1..dd] do Append(~vs,u);u:=u*Am;end for;
  ww:=Rank(Matrix(vs));
  matrixfile:=Sprintf("work/deg10/exclusions/krylov_N%o_chi%o_order%o.json",N,chiindex,Order(chi));
  Write(matrixfile,Sprintf("{\"N\":%o,\"chi_index\":%o,\"order\":%o,\"p\":%o,\"d\":%o,\"v\":%o,\"A\":%o}",N,chiindex,Order(chi),#F,dd,[Integers()!x:x in Eltseq(vm)],[Integers()!x:x in Eltseq(Am)]));
  printf "WINDING_KRYLOV_MATRIX N=%o chi_index=%o order=%o path=%o\n",N,chiindex,Order(chi),matrixfile;
  printf "WINDING_KRYLOV N=%o order=%o span_rank=%o dimension=%o verified=%o time=%o\n",N,Order(chi),ww,dd,ww eq dd,Cputime(tt);
  if ww ne dd then ww:=Dimension(WindingSubmodule(S));end if;
 else
  ww:=Dimension(WindingSubmodule(S));
 end if;
 printf "WINDING N=%o order=%o conductor=%o relative_dimension=%o winding_dimension=%o coefficient_degree=%o time=%o\n",N,Order(chi),Conductor(chi),dd,ww,Degree(BaseRing(S)),Cputime(tt);
 total+:=dd*EulerPhi(Order(chi));
 Append(~allchi,chi);if dd ne ww then Append(~bad,chi);end if;
end for;
I:=EulerPhi(N)*N*&*[Rationals()|1+1/p[1]:p in Factorization(N)]/2;
g:=Integers()!(1+I/12-&+[EulerPhi(d)*EulerPhi(N div d):d in Divisors(N)]/4);
assert total eq g;
admissible:=[a:a in [2..N-2]|GCD(a,N) eq 1 and &and[Evaluate(chi,a) eq 1:chi in bad]];
printf "WINDING_FINAL N=%o genus=%o total_dimension=%o bad_orders=%o admissible=%o\n",N,g,total,[Order(chi):chi in bad],admissible;
print "WINDING_CHARACTER_PROBE_DONE";
quit;
