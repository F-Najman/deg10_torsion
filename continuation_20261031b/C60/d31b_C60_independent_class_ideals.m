SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(6*10^9);Alarm(850);tt:=Cputime();AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");root:="work/deg10/continuation_20261031b/C60/";a:=eval Read(root cat "target30_proven_picard.m");X:=MDX1(30,GF(7));C:=Curve(X);F0,mp:=AlgorithmicFunctionField(FunctionField(C));F:=RationalExtensionRepresentation(F0);assert Genus(F)eq 9;k:=GF(7);PK:=PolynomialRing(k);B:=BaseField(F);
function Rat(v)return B!(PK!v[1])/B!(PK!v[2]);end function;
assert Coefficients(DefiningPolynomial(F))eq[Rat(v):v in a[1]];
function FF(v)return F![Rat(z):z in v];end function;
xy:=[F!mp(v):v in X`_coordinates_xy];assert xy eq[FF(v):v in a[7]];OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
function DecodePrime(v)
 O:=v[2]eq 1 select OF else OI;J:=ideal<O|[FF(z):z in v[3]]>;assert IsPrime(J);P:=Place(J);assert Degree(P)eq v[1];return P;
end function;
fb:=[DecodePrime([*a[3][i],a[2][i],a[4][i]*]):i in[1..#a[3]]];AB:=DecodePrime(a[12]);assert Degree(AB)eq 1;A:=Divisor(AB);ia:=Index(fb,AB);assert ia gt 0;Z:=FreeAbelianGroup(#a[6][1]);G,q:=quo<Z|[Z!v:v in a[6]]>;fc:=[q(Z!v):v in a[5]];assert #fb eq 164;
function FromKey(key)
 e:=key[1];bits:=key[2..3];pol:=PK!key[4..e+4];assert Degree(pol)eq e and IsIrreducible(pol);cs:=[PK!key[e+5..2*e+4],PK!key[2*e+5..3*e+4]];g:=bits[1]*xy[1]+bits[2]*xy[2];fs:=[Evaluate(pol,g),xy[1]-Evaluate(cs[1],g),xy[2]-Evaluate(cs[2],g)];J:=ideal<OF|fs>;ff:=Factorization(J);pp:=[Place(v[1]):v in ff|v[2]gt 0];if #pp eq 0 then J:=ideal<OI|fs>;ff:=Factorization(J);pp:=[Place(v[1]):v in ff|v[2]gt 0];end if;assert #pp eq 1;P:=pp[1];assert Degree(P)eq e and &and[f eq 0 or Valuation(f,P)gt 0:f in fs];vv:=[Evaluate(z,P):z in xy];gg:=bits[1]*vv[1]+bits[2]*vv[2];assert MinimalPolynomial(gg,k)eq pol and &and[vv[j]eq Evaluate(cs[j],gg):j in[1,2]];return P;
end function;
function SparseDiv(v)
 D:=DivisorGroup(F)!0;for z in v do assert z[1]in[1..#fb];D+:=z[2]*Divisor(fb[z[1]]);end for;return D;
end function;
function SparseClass(v)
 return &+[G|z[2]*fc[z[1]]:z in v];
end function;
out:=[* *];counts:=[];nt:=0;
for shard in[0..2]do packet:=eval Read(root cat Sprintf("three_maps_s%o_v2.m",shard));W:=eval Read(root cat Sprintf("class_witnesses_s%o.m",shard));wi:=0;direct:=0;rs:=[* *];
 for ix in[1..#packet[2]]do im:=packet[2][ix];P:=FromKey(im[1]);e:=Degree(P);assert im[3]eq e;ii:=e le Maximum(a[3]) select Index(fb,P)else 0;target:=q(Z!im[2]);
  if ii gt 0 then assert target eq fc[ii]-e*fc[ia];direct+:=1;used:=0;
  else wi+:=1;assert wi le #W;v:=W[wi];assert #v eq 7 and v[2]eq 0 and q(Z!v[3])eq target;E:=SparseDiv(v[4]);R:=SparseDiv(v[5]);AA:=SparseDiv(v[6]);assert E eq DivisorGroup(F)!E and Degree(E)eq 9 and Degree(R)eq 9 and AA eq A;assert &and[z[2]ge 0:z in v[4]]and &and[z[2]ge 0:z in v[5]];b:=FF(v[7]);assert b ne 0;IE,JE:=Ideals(E);IR,JR:=Ideals(R);ID,JD:=Ideals(Divisor(P)-e*A);assert IE/(IR*(b*OF))eq ID and JE/(JR*(b*OI))eq JD;assert SparseClass(v[4])-SparseClass(v[5])eq target;used:=wi;
  end if;Append(~rs,[*ix,e,ii,used,Eltseq(target)*]);nt+:=1;if ix mod 100 eq 0 then print "D31B_C60_NATIVE_IDEAL_REPLAY_PROGRESS",shard,ix,wi,direct,Cputime(tt);end if;
 end for;assert wi eq #W and #packet[2]eq 837;Append(~out,rs);Append(~counts,[shard,#rs,wi,direct]);print "D31B_C60_NATIVE_IDEAL_REPLAY_SHARD_COMPLETE",shard,#rs,wi,direct,Cputime(tt);
end for;
assert nt eq 2511;Write(root cat "independent_native_class_ideals.m",Sprint([*out,counts*]):Overwrite:=true);print "D31B_C60_NATIVE_IDEAL_REPLAY_COMPLETE",nt,counts,Cputime(tt);quit;
