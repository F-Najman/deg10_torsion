SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(4*10^9);tt:=Cputime();AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
Fp:=GF(3);PR:=PolynomialRing(Fp);
for D in DEGREES do
 K:=GF(3^D);PK:=PolynomialRing(K);seenj:={K|};orb56:=AssociativeArray();orb70:=AssociativeArray();points56:=AssociativeArray();points70:=AssociativeArray();jcnt:=0;curvecnt:=0;qual56:=0;qual70:=0;
 function Ky(v)
  seq:=[];for i in[1..D] do Append(~seq,&cat[[Integers()!a:a in Eltseq(z,Fp)]:z in v]);v:=[z^3:z in v];end for;return Sort(seq)[1];
 end function;
 function Deg(v)
  for e in Divisors(D) do if &and[z^(3^e) eq z:z in v]then return e;end if;end for;error "fielddegree";
 end function;
 for j in K do
  if j in seenj then continue;end if;for i in[0..D-1]do Include(~seenj,j^(3^i));end for;jcnt+:=1;es:=Twists(EllipticCurveWithjInvariant(j));
  for E in es do
   curvecnt+:=1;card:=#E;if card mod 56 ne 0 and card mod 70 ne 0 then continue;end if;G,gm:=AbelianGroup(E);ivs:=Invariants(G);assert #ivs in[1,2];
   for N in[56,70]do
    if Exponent(G) mod N ne 0 then continue;end if;if N eq 56 then qual56+:=1;else qual70+:=1;end if;
    ms:=[GCD(v,N):v in ivs];gg:=[(ivs[i] div ms[i])*G.i:i in[1..#ivs]];tors:=#ivs eq 1 select [a*gg[1]:a in[0..ms[1]-1]] else [a*gg[1]+b*gg[2]:a in[0..ms[1]-1],b in[0..ms[2]-1]];done:={G|};
    for g in tors do
     if N*g ne G!0 or exists{p:p in PrimeDivisors(N)|(N div p)*g eq G!0}or g in done then continue;end if;done join:={u*g:u in[1..N-1]|GCD(u,N)eq 1};P:=gm(g);assert Order(P)eq N;
     ks:={Ky(TateNormalForm_bc(E,u*P)):u in[1..N div 2]|GCD(u,N)eq 1};ok:=Sort(Setseq(ks))[1];v:=[Evaluate(PR!ok[1..D],K.1),Evaluate(PR!ok[D+1..2*D],K.1)];e:=Deg(v);
     if N eq 56 then if not IsDefined(orb56,ok)then orb56[ok]:=e;end if;for k in ks do points56[k]:=e;end for;
     else if not IsDefined(orb70,ok)then orb70[ok]:=e;end if;for k in ks do points70[k]:=e;end for;end if;
    end for;
   end for;
  end for;
  if jcnt mod 500 eq 0 then print "D31A_P3_SOURCE_PROGRESS",D,jcnt,curvecnt,qual56,qual70,#Keys(points56),#Keys(points70),Cputime(tt);end if;
 end for;assert #seenj eq #K;
 for N in[56,70]do
  oo:=N eq 56 select orb56 else orb70;pp:=N eq 56 select points56 else points70;recs:=[[*oo[v],v*]:v in Sort(Setseq(Keys(oo)))];cts:=[#[v:v in Keys(pp)|pp[v]eq e]:e in[1..D]];
  Write(Sprintf("work/deg10/continuation_20261031a/C%o/source_p3_d%o.m",N,D),Sprint([*D,[Integers()!a:a in Coefficients(DefiningPolynomial(K,Fp))],true,jcnt,curvecnt,cts,recs*]):Overwrite:=true);print "D31A_P3_SOURCE_DEGREE_COMPLETE",N,D,cts,#recs,Cputime(tt);
 end for;
end for;print "D31A_P3_ALL_SOURCE_CENSUS_COMPLETE",DEGREES,Cputime(tt);quit;
