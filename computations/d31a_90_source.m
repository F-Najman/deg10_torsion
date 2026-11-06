SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(4*10^9);Alarm(240);tt:=Cputime();AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");root:="work/deg10/continuation_20261031a/C90/";data:=eval Read("work/deg10/continuation_20261031a/C90/inputs/all_curves.m");k:=GF(2);R:=PolynomialRing(k);fields:=[* *];coverage:=[* *];records:=[* *];aa:=[a:a in [1..22]|GCD(a,45) eq 1];assert #aa eq EulerPhi(45) div 2;
function Key(v,D) keys:=[];for i in [1..D] do Append(~keys,&cat[[Integers()!x:x in Eltseq(z,GF(2))]:z in v]);v:=[z^2:z in v];end for;return Sort(keys)[1];end function;
for dat in data do
 if dat[1] ne 10 then continue;end if;
 d,pol,curves:=Explode(dat);K<z>:=ext<k|R!pol>;Append(~fields,[*d,pol*]);points:={};orbs:=AssociativeArray();groups:=[* *];qual:=0;cycs:=0;
 for er in curves do
  if er[3] mod 45 ne 0 then continue;end if;E:=EllipticCurve([Evaluate(R!v,z):v in er[2]]);assert #E eq er[3];G,gm:=AbelianGroup(E);Append(~groups,Invariants(G));if Exponent(G) mod 45 ne 0 then continue;end if;qual+:=1;seen:={G|};
  for g in G do
   if 45*g ne G!0 or exists{p:p in PrimeDivisors(45)|(45 div p)*g eq G!0} or g in seen then continue;end if;P:=gm(g);assert 45*P eq E!0 and &and[(45 div p)*P ne E!0:p in PrimeDivisors(45)];cycs+:=1;seen join:={a*g:a in [1..44]|GCD(a,45) eq 1};assert TateNormalForm_bc(E,P) eq TateNormalForm_bc(E,-P);ks:={Key(TateNormalForm_bc(E,a*P),d):a in aa};points join:=ks;ky:=Sort(Setseq(ks))[1];
   if not IsDefined(orbs,ky) then orbs[ky]:=true;end if;
  end for;
 end for;
 keys:=Sort(Setseq(Keys(orbs)));for ky in keys do b:=Evaluate(R!ky[1..d],z);c:=Evaluate(R!ky[d+1..2*d],z);EC:=EllipticCurve([K|1-c,-b,-b,0,0]);P:=EC![0,0];assert 45*P eq EC!0 and &and[(45 div p)*P ne EC!0:p in PrimeDivisors(45)];assert &and[exists{w:w in [b,c]|w^(2^e) ne w}:e in Divisors(d)|e lt d];Append(~records,[*d,ky*]);end for;
  Append(~coverage,[*d,#curves,qual,groups,cycs,#points,#keys*]);print "C90_DYADIC_SOURCE_DEGREE_COMPLETE",d,#curves,qual,#points,#keys,groups,Cputime(tt);
end for;
Write(root cat "source45.m",Sprint([*fields,coverage,records*]):Overwrite:=true);print "C90_DYADIC_SOURCE_CENSUS_COMPLETE",#records,Cputime(tt);quit;
