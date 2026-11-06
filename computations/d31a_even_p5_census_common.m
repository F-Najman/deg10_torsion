SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(6*10^9);tt:=Cputime();AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");Fp:=GF(5);PR:=PolynomialRing(Fp);LEVELS:=[48,54,66];
for D in DEGREES do
 K:=GF(5^D);seenj:={K|};orbs:=AssociativeArray();points:=AssociativeArray();jcnt:=0;curves:=0;selected:=0;complete:=true;
 function Ky(v)
  seq:=[];for i in[1..D]do Append(~seq,&cat[[Integers()!a:a in Eltseq(z,Fp)]:z in v]);v:=[z^5:z in v];end for;return Sort(seq)[1];
 end function;
 function Deg(v)
  for e in Divisors(D)do if &and[z^(5^e)eq z:z in v]then return e;end if;end for;error "fielddegree";
 end function;
 for j in K do
  if j in seenj then continue;end if;if LIMIT gt 0 and jcnt ge LIMIT then complete:=false;break;end if;for i in[0..D-1]do Include(~seenj,j^(5^i));end for;jcnt+:=1;if(jcnt-1)mod NSHARDS ne SHARD then continue;end if;selected+:=1;
  for E in Twists(EllipticCurveWithjInvariant(j))do
   curves+:=1;card:=#E;eligible:=[N:N in LEVELS|card mod N eq 0];if #eligible eq 0 then continue;end if;G,gm:=AbelianGroup(E);ivs:=Invariants(G);assert #ivs in[1,2];
   for N in eligible do
    if Exponent(G)mod N ne 0 then continue;end if;ms:=[GCD(v,N):v in ivs];gg:=[(ivs[i]div ms[i])*G.i:i in[1..#ivs]];tors:=#ivs eq 1 select[a*gg[1]:a in[0..ms[1]-1]]else[a*gg[1]+b*gg[2]:a in[0..ms[1]-1],b in[0..ms[2]-1]];done:={G|};
    for g in tors do
     if N*g ne G!0 or exists{p:p in PrimeDivisors(N)|(N div p)*g eq G!0}or g in done then continue;end if;done join:={u*g:u in[1..N-1]|GCD(u,N)eq 1};P:=gm(g);assert Order(P)eq N;ks:={Ky(TateNormalForm_bc(E,u*P)):u in[1..N div 2]|GCD(u,N)eq 1};ok:=Sort(Setseq(ks))[1];v:=[Evaluate(PR!ok[1..D],K.1),Evaluate(PR!ok[D+1..2*D],K.1)];e:=Deg(v);orbs[[N]cat ok]:=e;for ky in ks do points[[N]cat ky]:=e;end for;
    end for;
   end for;
  end for;
  if selected mod 1000 eq 0 then print "D31A_P5_SOURCE_PROGRESS",D,SHARD,NSHARDS,jcnt,selected,curves,#Keys(orbs),#Keys(points),Cputime(tt);end if;
 end for;
 if complete then assert #seenj eq #K;end if;
 for N in LEVELS do
  okeys:=Sort([v:v in Keys(orbs)|v[1]eq N]);recs:=[[*orbs[v],v[2..#v]*]:v in okeys];cts:=[#[v:v in Keys(points)|v[1]eq N and points[v]eq e]:e in[1..D]];
  name:=LIMIT gt 0 select "source_p5_profile1000.m"else Sprintf("source_p5_d%o_s%o.m",D,SHARD);Write(Sprintf("work/deg10/continuation_20261031a/C%o/",N)cat name,Sprint([*D,[Integers()!a:a in Coefficients(DefiningPolynomial(K,Fp))],complete,SHARD,NSHARDS,jcnt,selected,curves,cts,recs*]):Overwrite:=true);print "D31A_P5_SOURCE_SHARD_DONE",N,D,SHARD,complete,cts,#recs,Cputime(tt);
 end for;
end for;print "D31A_P5_CENSUS_RUN_FINISHED",DEGREES,SHARD,NSHARDS,LIMIT,Cputime(tt);quit;
