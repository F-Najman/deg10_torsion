SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(4*10^9);tt:=Cputime();AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");Attach("work/deg10/continuation_20261030c/C91_CM_control/inputs/fast_hecke.m");
root:="work/deg10/continuation_20261030c/C91_CM_control/";load "work/deg10/continuation_20261030c/C91_CM_control/inputs/records91_8.m";k:=GF(2);PK:=PolynomialRing(k);K<z>:=ext<k|PK!FieldPolynomial>;
function Key(v)
 L:=Parent(v[1]);d:=LCM([Degree(MinimalPolynomial(a,GF(2))):a in v]);found:=false;chosen:=[Integers()|];
 for bits in [[1,0],[0,1],[1,1]] do g:=&+[L!bits[i]*v[i]:i in[1,2]];po:=MinimalPolynomial(g,GF(2));if Degree(po) eq d then found:=true;chosen:=bits;break;end if;end for;assert found;
 M:=Matrix(GF(2),[Eltseq(g^i,GF(2)):i in[0..d-1]]);assert Rank(M) eq d;co:=[Eltseq(Solution(M,Vector(GF(2),Eltseq(a,GF(2))))):a in v];return [*d,chosen,[Integers()!a:a in Coefficients(po)],[[Integers()!a:a in row]:row in co]*];
end function;
function PointKey(E,P)return Key(TateNormalForm_xy(E,P));end function;
ss:=[* *];seen:={};
for rec in Records do E:=EllipticCurve([Evaluate(PK!v,z):v in rec[2]]);P:=E![Evaluate(PK!v,z):v in rec[3]];assert jInvariant(E) eq 0 and #E eq 273 and Order(P) eq 91;
 for u in [u:u in[1..45]|GCD(u,91) eq 1] do Q:=u*P;ky:=PointKey(E,Q);assert ky[1] eq 8 and PointKey(E,16*Q) eq ky;if not Sprint(ky) in seen then Include(~seen,Sprint(ky));Append(~ss,[*E,Q,ky*]);end if;end for;
end for;assert #ss eq 3;
keys:=[r[3]:r in ss];delta:=[Index(keys,PointKey(r[1],9*r[2])):r in ss];assert Sort(delta) eq[1..3] and &and[delta[i] ne i:i in[1..3]];
out:=[* *];
for l in PRIMES do
 for si in[1..3] do E,Q,ky:=Explode(ss[si]);fast:=MDIsogeniesFast(E,Q,l);assert &+[v[3]:v in fast] eq l+1;trs:=[* *];
  for tr in fast do ET,QT,e:=Explode(tr);assert jInvariant(ET) eq 0 and Order(QT) eq 91;tk:=PointKey(ET,QT);assert (8*e) mod tk[1] eq 0;mult:=8*e div tk[1];Append(~trs,[*e,mult,tk,PointKey(ET,9*QT),PointKey(ET,81*QT)*]);end for;
  for a in[9,81] do
   coeff:=AssociativeArray();dgs:=AssociativeArray();
   for tr in trs do
    for term in [[*tr[a eq 9 select 4 else 5],tr[2]*],[*tr[3],-tr[2]*]] do sk:=Sprint(term[1]);if not IsDefined(coeff,sk) then coeff[sk]:=0;dgs[sk]:=term[1][1];end if;coeff[sk]+:=term[2];end for;
   end for;
   for term in [[*PointKey(E,a*Q),-1*],[*ky,1*],[*PointKey(E,a*l*Q),-l*],[*PointKey(E,l*Q),l*]] do sk:=Sprint(term[1]);if not IsDefined(coeff,sk) then coeff[sk]:=0;dgs[sk]:=term[1][1];end if;coeff[sk]+:=term[2];end for;
   nz:=[[*eval sk,coeff[sk]*]:sk in Keys(coeff)|coeff[sk] ne 0];pos:=&+[Integers()|dgs[sk]*Maximum(0,coeff[sk]):sk in Keys(coeff)];neg:=&+[Integers()|dgs[sk]*Maximum(0,-coeff[sk]):sk in Keys(coeff)];assert pos eq neg;
   Append(~out,[*l,si,a,ky,trs,nz,pos,neg*]);print "C91_CM_MARKED_HECKE_CASE",l,si,a,[r[1]:r in trs],[r[3][1]:r in trs],[r[2]:r in trs],pos,neg,Cputime(tt);
  end for;
 end for;
end for;
Write(root cat Sprintf("Hecke_CM_%o.m",PRIMES[1]),Sprint([*keys,delta,out*]):Overwrite:=true);print "C91_CM_MARKED_HECKE_COMPLETE",PRIMES,#out,Cputime(tt);quit;
