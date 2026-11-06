SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);
AttachSpec("work/degree9/shared/mdmagma/v2/mdmagma.spec");Attach("work/degree9/shared/fast_hecke.m");
tt:=Cputime();
function XYKey(E,P)
 xx,yy:=Explode(TateNormalForm_xy(E,P));
 d:=LCM(Degree(MinimalPolynomial(xx,GF(2))),Degree(MinimalPolynomial(yy,GF(2))));
 found:=false;
 for bits in [[1,0],[0,1],[1,1]] do
  z:=bits[1]*xx+bits[2]*yy;mp:=MinimalPolynomial(z,GF(2));
  if Degree(mp) eq d then bb:=bits;found:=true;break;end if;
 end for;assert found;
 B:=Matrix(GF(2),[Eltseq(z^i,GF(2)):i in [0..d-1]]);assert Rank(B) eq d;
 cs:=[Eltseq(Solution(B,Vector(GF(2),Eltseq(v,GF(2))))):v in [xx,yy]];
 return [d] cat bb cat [Integers()!c:c in Eltseq(mp)] cat [Integers()!c:c in &cat cs];
end function;
function Reconstruct(r)
 d,i,fp,ecs,pcs:=Explode(r);R<z>:=PolynomialRing(GF(2));F<a>:=ext<GF(2)|R!fp>;
 assert Degree(F) eq d;
 fromvec:=func<v|&+[F!v[k]*a^(k-1):k in [1..#v]]>;
 E:=EllipticCurve([fromvec(v):v in ecs]);P:=E![fromvec(v):v in pcs];
 assert Order(P) eq N;return E,P;
end function;
for rec in Inputs do
 d:=rec[1];i:=rec[2];E,P:=Reconstruct(rec);inputkey:=XYKey(E,P);assert inputkey[1] eq d;
 phs:=MDIsogenies(E,5);assert #phs eq 6;
 native:=AssociativeArray();fast:=AssociativeArray();native_a:=AssociativeArray();fast_a:=AssociativeArray();
 for ph in phs do
  EC:=Codomain(ph);PC:=ph(P);assert Order(PC) eq N;k:=XYKey(EC,PC);
  if not IsDefined(native,k) then native[k]:=0;end if;native[k]+:=1;
  if a ne 1 then ka:=XYKey(EC,a*PC);if not IsDefined(native_a,ka) then native_a[ka]:=0;end if;native_a[ka]+:=1;end if;
 end for;
 ins:=MDIsogeniesFast(E,P,5);assert &+[z[3]:z in ins] eq 6;
 for z in ins do
  EC,PC,m:=Explode(z);assert Order(PC) eq N;k:=XYKey(EC,PC);
  if not IsDefined(fast,k) then fast[k]:=0;end if;fast[k]+:=m;
  if a ne 1 then ka:=XYKey(EC,a*PC);if not IsDefined(fast_a,ka) then fast_a[ka]:=0;end if;fast_a[ka]+:=m;end if;
 end for;
 assert Keys(native) eq Keys(fast);
 for k in Keys(native) do assert native[k] eq fast[k] and (d*native[k]) mod k[1] eq 0;end for;
 if a ne 1 then assert Keys(native_a) eq Keys(fast_a);for k in Keys(native_a) do assert native_a[k] eq fast_a[k] and (d*native_a[k]) mod k[1] eq 0;end for;end if;
 printf "DIRECT69_INDEPENDENT_HECKE_CONTROL N=%o degree=%o index=%o native_maps=6 fast_orbits=%o weighted_image_degrees=%o observations=%o weights=%o alpha=%o agrees=true time=%o\n",N,d,i,#ins,[k[1]:k in Keys(native)],[native[k]:k in Keys(native)],[d*native[k] div k[1]:k in Keys(native)],a,Cputime(tt);
end for;
printf "DIRECT69_NATIVE_HECKE_CONTROL_COMPLETE N=%o inputs=%o time=%o\n",N,#Inputs,Cputime(tt);quit;
