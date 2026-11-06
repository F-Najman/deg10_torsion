SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(5*10^9);Alarm(500);tt:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
n:=33;root:="work/deg10/continuation_20261031b/M3N33/";data:=eval Read(root cat "inputs/all_curves.m");dat:=[v:v in data|v[1] eq 10][1];k:=GF(2);PP:=PolynomialRing(k);K<z>:=ext<k|PP!dat[2]>;zet:=RootOfUnity(3,K);assert zet^4 eq zet;
function Key(v)
 keys:=[];for i in[1..5] do Append(~keys,&cat[[Integers()!a:a in Eltseq(w,k)]:w in v]);v:=[w^4:w in v];end for;return Sort(keys)[1];end function;
function Normal3(E,P,Q)
 a1,a2,a3,a4,a6:=Explode(aInvariants(E));r:=P[1];tt:=P[2];sl:=(3*r^2+2*a2*r+a4-a1*tt)/(2*tt+a1*r+a3);aa1:=a1+2*sl;aa3:=a3+r*a1+2*tt;xx:=Q[1]-r;yy:=Q[2]-tt-sl*xx;sc:=yy/xx^2;t:=-sc^2*xx;q:=(sc*aa1-t-2)/(t-1);assert sc^3*aa3 eq q*t^2+(1-q)*t;return [q,t];end function;
allkeys:={};cyclics:={};representatives:=[* *];coverage:=[* *];seen:={};aa:=[a:a in[1..n-1]|GCD(a,n) eq 1];
for ix in[1..#dat[3]] do
 er:=dat[3][ix];if er[3] mod (3*n) ne 0 then continue;end if;E:=EllipticCurve([Evaluate(PP!v,z):v in er[2]]);assert #E eq er[3];G,gm:=AbelianGroup(E);inv:=Invariants(G);if #inv ne 2 or inv[1] mod 3 ne 0 or inv[2] mod n ne 0 then continue;end if;
 qs:=[gm(g):g in G|Order(g) eq n];ps:=[gm(g):g in G|Order(g) eq 3];localkeys:={};qkeys:={};
 for Q in qs do
  ky:=Key(TateNormalForm_bc(E,Q));Include(~qkeys,ky);Include(~cyclics,ky);
  orbkeys:={Key(TateNormalForm_bc(E,a*Q)):a in aa};orb:=Sort(Setseq(orbkeys))[1];
  if orb notin seen then Include(~seen,orb);Append(~representatives,[*ix,inv,er,orb,Setseq(orbkeys),#BaseRing(E),jInvariant(E) eq 0 select 1 else 0*]);end if;
  for P in ps do
   if WeilPairing(P,(n div 3)*Q,3) ne zet then continue;end if;
   ky3:=Key(Normal3(E,P,Q));Include(~allkeys,ky3);Include(~localkeys,ky3);
  end for;
 end for;
 Append(~coverage,[*ix,inv,#localkeys,#qkeys,er[3]*]);print "RELATIVE5_PHYSICAL_CURVE",n,ix,inv,#localkeys,#qkeys,#seen,Cputime(tt);
end for;
Write(root cat "physical_source5.m",Sprint([*dat[2],Eltseq(zet,k),coverage,Sort(Setseq(allkeys)),Sort(Setseq(cyclics)),representatives*]):Overwrite:=true);print "D31B_PHYSICAL_SOURCE5_COMPLETE",n,#allkeys,#cyclics,#representatives,Cputime(tt);quit;
