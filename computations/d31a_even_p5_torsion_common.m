SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(6*10^9);tt:=Cputime();AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");Fp:=GF(5);PR:=PolynomialRing(Fp);LEVELS:=[48,54,66];
for D in DEGREES do
 K:=GF(5^D);PX<X>:=PolynomialRing(K);seenj:={K|};orbs:=AssociativeArray();points:=AssociativeArray();jcnt:=0;curves:=0;selected:=0;complete:=true;good2:=0;good3:=0;
 function Code(v)return[Integers()!a:a in Eltseq(v,Fp)];end function;
 function Ky(v)seq:=[];for i in[1..D]do Append(~seq,&cat[Code(z):z in v]);v:=[z^5:z in v];end for;return Sort(seq)[1];end function;
 function Deg(v)for e in Divisors(D)do if &and[z^(5^e)eq z:z in v]then return e;end if;end for;error "fielddegree";end function;
 function Subkey(Q,m)return Sort([Code((u*Q)[1]):u in[1..m div 2]|GCD(u,m)eq 1]);end function;
 function Unique(pp,m)ss:={};out:=[* *];for Q in pp do ky:=Subkey(Q,m);if ky in ss then continue;end if;Include(~ss,ky);Append(~out,Q);end for;return out;end function;
 function Halves(E,R)
  b2,b4,b6,b8:=Explode(bInvariants(E));ai:=aInvariants(E);cubic:=X^3+ai[2]*X^2+ai[4]*X+ai[5];pol:=X^4-b4*X^2-2*b6*X-b8-R[1]*(4*X^3+b2*X^2+2*b4*X+b6);out:=[* *];
  for rr in Roots(pol)do xx:=rr[1];vv:=Evaluate(cubic,xx);if not IsSquare(vv)then continue;end if;Q:=E![xx,Sqrt(vv)];if 2*Q eq -R then Q:=-Q;end if;assert 2*Q eq R;Append(~out,Q);end for;return out;
 end function;
 function Thirds(E,R)
  b2,b4,b6,b8:=Explode(bInvariants(E));ai:=aInvariants(E);cubic:=X^3+ai[2]*X^2+ai[4]*X+ai[5];p3:=DivisionPolynomial(E,3);p4b:=2*X^6+b2*X^5+5*b4*X^4+10*b6*X^3+10*b8*X^2+(b2*b8-b4*b6)*X+b4*b8-b6^2;num:=X*p3^2-4*cubic*p4b;pol:=num-R[1]*p3^2;assert Degree(pol)eq 9;out:=[* *];
  for rr in Roots(pol)do xx:=rr[1];vv:=Evaluate(cubic,xx);if not IsSquare(vv)then continue;end if;Q:=E![xx,Sqrt(vv)];if 3*Q eq -R then Q:=-Q;end if;assert 3*Q eq R;Append(~out,Q);end for;return out;
 end function;
 for j in K do
  if j in seenj then continue;end if;if LIMIT gt 0 and selected ge LIMIT then complete:=false;break;end if;for i in[0..D-1]do Include(~seenj,j^(5^i));end for;jcnt+:=1;if(jcnt-1)mod NSHARDS ne SHARD then continue;end if;if SKIPPROPER and Degree(MinimalPolynomial(j,Fp))ne D then continue;end if;selected+:=1;es:=Twists(EllipticCurveWithjInvariant(j));
  for EE in es do
   curves+:=1;ai:=aInvariants(EE);E:=EllipticCurve([K|0,ai[2]+ai[1]^2/4,0,ai[4]+ai[1]*ai[3]/2,ai[5]+ai[3]^2/4]);ai:=aInvariants(E);cubic:=X^3+ai[2]*X^2+ai[4]*X+ai[5];twos:=[*E![r[1],0]:r in Roots(cubic)*];if #twos eq 0 then if #es eq 2 then break;else continue;end if;end if;good2+:=1;
   threes:=[* *];for rr in Roots(DivisionPolynomial(E,3))do xx:=rr[1];vv:=Evaluate(cubic,xx);if IsSquare(vv)then Q:=E![xx,Sqrt(vv)];assert 3*Q eq E!0 and Q ne E!0;Append(~threes,Q);end if;end for;if #threes eq 0 then continue;end if;good3+:=1;
   fours:=[* *];for R in twos do fours cat:=Halves(E,R);end for;fours:=Unique(fours,4);eights:=[* *];for R in fours do eights cat:=Halves(E,R);end for;eights:=Unique(eights,8);sixteens:=[* *];for R in eights do sixteens cat:=Halves(E,R);end for;sixteens:=Unique(sixteens,16);
   for Q16 in sixteens do for Q3 in threes do Q:=Q16+Q3;assert 48*Q eq E!0 and 24*Q ne E!0 and 16*Q ne E!0;N:=48;ks:={Ky(TateNormalForm_bc(E,u*Q)):u in[1..24]|GCD(u,48)eq 1};ok:=Sort(Setseq(ks))[1];v:=[Evaluate(PR!ok[1..D],K.1),Evaluate(PR!ok[D+1..2*D],K.1)];orbs[[N]cat ok]:=Deg(v);for ky in ks do points[[N]cat ky]:=Deg(v);end for;end for;end for;
   nines:=[* *];for R in threes do nines cat:=Thirds(E,R);end for;nines:=Unique(nines,9);twentysevens:=[* *];for R in nines do twentysevens cat:=Thirds(E,R);end for;twentysevens:=Unique(twentysevens,27);
   for Q27 in twentysevens do for Q2 in twos do Q:=Q27+Q2;assert 54*Q eq E!0 and 27*Q ne E!0 and 18*Q ne E!0;N:=54;ks:={Ky(TateNormalForm_bc(E,u*Q)):u in[1..27]|GCD(u,54)eq 1};ok:=Sort(Setseq(ks))[1];v:=[Evaluate(PR!ok[1..D],K.1),Evaluate(PR!ok[D+1..2*D],K.1)];orbs[[N]cat ok]:=Deg(v);for ky in ks do points[[N]cat ky]:=Deg(v);end for;end for;end for;
   elevens:=[* *];for rr in Roots(DivisionPolynomial(E,11))do xx:=rr[1];vv:=Evaluate(cubic,xx);if IsSquare(vv)then Q:=E![xx,Sqrt(vv)];assert 11*Q eq E!0 and Q ne E!0;Append(~elevens,Q);end if;end for;elevens:=Unique(elevens,11);
   for Q11 in elevens do for Q3 in threes do for Q2 in twos do Q:=Q11+Q3+Q2;assert 66*Q eq E!0 and 33*Q ne E!0 and 22*Q ne E!0 and 6*Q ne E!0;N:=66;ks:={Ky(TateNormalForm_bc(E,u*Q)):u in[1..33]|GCD(u,66)eq 1};ok:=Sort(Setseq(ks))[1];v:=[Evaluate(PR!ok[1..D],K.1),Evaluate(PR!ok[D+1..2*D],K.1)];orbs[[N]cat ok]:=Deg(v);for ky in ks do points[[N]cat ky]:=Deg(v);end for;end for;end for;end for;
  end for;
  if selected mod 1000 eq 0 then print "D31A_P5_TORSION_PROGRESS",D,SHARD,NSHARDS,jcnt,selected,curves,good2,good3,#Keys(orbs),#Keys(points),Cputime(tt);end if;
 end for;
 if complete then assert #seenj eq #K;end if;
 for N in LEVELS do
  okeys:=Sort([v:v in Keys(orbs)|v[1]eq N]);recs:=[[*orbs[v],v[2..#v]*]:v in okeys];cts:=[#[v:v in Keys(points)|v[1]eq N and points[v]eq e]:e in[1..D]];name:=LIMIT gt 0 select "source_torsion_profile.m"else Sprintf("source_torsion_p5_d%o_s%o.m",D,SHARD);Write(Sprintf("work/deg10/continuation_20261031a/C%o/",N)cat name,Sprint([*D,[Integers()!a:a in Coefficients(DefiningPolynomial(K,Fp))],complete,SHARD,NSHARDS,jcnt,selected,curves,cts,recs*]):Overwrite:=true);print "D31A_P5_TORSION_SOURCE_SHARD_DONE",N,D,SHARD,complete,cts,#recs,Cputime(tt);
 end for;
end for;print "D31A_P5_TORSION_RUN_FINISHED",DEGREES,SHARD,NSHARDS,LIMIT,Cputime(tt);quit;
