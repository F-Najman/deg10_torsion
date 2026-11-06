SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(4*10^9);tt:=Cputime();AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");root:="work/deg10/continuation_20261030b/C81/";data:=eval Read(root cat "inputs/all_curves.m");k:=GF(2);R:=PolynomialRing(k);fields:=[* *];coverage:=[* *];records:=[* *];aa:=[a:a in [1..40]|GCD(a,81) eq 1];assert #aa eq 27;
function Key(v,D) keys:=[];for i in [1..D] do Append(~keys,&cat[[Integers()!x:x in Eltseq(z,GF(2))]:z in v]);v:=[z^2:z in v];end for;return Sort(keys)[1];end function;
for dat in data do
 d,pol,curves:=Explode(dat);K<z>:=ext<k|R!pol>;Append(~fields,[*d,pol*]);points:={};orbs:=AssociativeArray();groups:=[* *];qual:=0;cycs:=0;
 for er in curves do
  if er[3] mod 81 ne 0 then continue;end if;E:=EllipticCurve([Evaluate(R!v,z):v in er[2]]);assert #E eq er[3];G,gm:=AbelianGroup(E);Append(~groups,Invariants(G));if Exponent(G) mod 81 ne 0 then continue;end if;qual+:=1;seen:={G|};
  for g in G do
   if 81*g ne G!0 or 27*g eq G!0 or g in seen then continue;end if;P:=gm(g);assert 81*P eq E!0 and 27*P ne E!0;cycs+:=1;seen join:={a*g:a in [1..80]|GCD(a,81) eq 1};assert TateNormalForm_bc(E,P) eq TateNormalForm_bc(E,-P);ks:={Key(TateNormalForm_bc(E,a*P),d):a in aa};points join:=ks;ky:=Sort(Setseq(ks))[1];
   if not IsDefined(orbs,ky) then orbs[ky]:=true;end if;
  end for;
 end for;
 keys:=Sort(Setseq(Keys(orbs)));for ky in keys do b:=Evaluate(R!ky[1..d],z);c:=Evaluate(R!ky[d+1..2*d],z);EC:=EllipticCurve([K|1-c,-b,-b,0,0]);P:=EC![0,0];assert 81*P eq EC!0 and 27*P ne EC!0;assert &and[exists{w:w in [b,c]|w^(2^e) ne w}:e in Divisors(d)|e lt d];Append(~records,[*d,ky*]);end for;
  Append(~coverage,[*d,#curves,qual,groups,cycs,#points,#keys*]);print "C81_DYADIC_SOURCE_DEGREE_COMPLETE",d,#curves,qual,#points,#keys,groups,Cputime(tt);
end for;
assert &and[r[6] eq 0:r in coverage|r[1] le 8];Write(root cat "source81.m",Sprint([*fields,coverage,records*]):Overwrite:=true);print "C81_DYADIC_SOURCE_CENSUS_COMPLETE",#records,Cputime(tt);quit;
