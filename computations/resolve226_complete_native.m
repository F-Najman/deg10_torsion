SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(3*10^9);tt:=Cputime();QQ:=Rationals();PQ<x>:=PolynomialRing(QQ);a:=eval Read("work/deg10/C226_resolution_20261030/complete_resolution_inputs.m");source:=eval Read("work/deg10/continuation_20261029b/C226/filter26.m");out:=[* *];
function Rat(v) return QQ!v[1]/QQ!v[2];end function;
function Ky(v,d)
 ims:=[];for k in [1..d] do Append(~ims,&cat[[Integers()!a:a in Eltseq(z,GF(3))]:z in v]);v:=[z^3:z in v];end for;return Sort(ims)[1];
end function;
for row in a do
 idx:=row[1];pt:=row[2];f:=PQ![Rat(v):v in pt[1]];assert Degree(f) eq 10 and IsIrreducible(f);K<alpha>:=NumberField(f);b:=K![Rat(v):v in pt[2]];c:=K![Rat(v):v in pt[3]];E:=EllipticCurve([K|0,c,0,b*(1-b-c),0]);P:=E![0,0];Q:=E![b,b];assert Discriminant(E) ne 0 and P ne E!0 and 2*P eq E!0 and 26*Q eq E!0 and 13*Q ne E!0 and 2*Q ne E!0 and P ne 13*Q;
 red:=row[3];p:=red[1];rr:=GF(p)!red[2];Rp:=PolynomialRing(GF(p));ff:=Rp!f;assert Degree(GCD(ff,Derivative(ff))) eq 0 and Evaluate(ff,rr) eq 0;bb:=Evaluate(Rp![Rat(v):v in pt[2]],rr);cc:=Evaluate(Rp![Rat(v):v in pt[3]],rr);Ep:=EllipticCurve([GF(p)|0,cc,0,bb*(1-bb-cc),0]);assert Discriminant(Ep) ne 0 and #Ep eq 52 and [Integers()!bb,Integers()!cc] eq red[3..4];
 d:=idx eq 30 select 5 else 10;R3:=PolynomialRing(GF(3));cov:=source[4][d-4];assert cov[1] eq d;L<z>:=ext<GF(3)|R3!cov[2]>;PL:=PolynomialRing(L);roots:=Roots(PL!f);assert #roots eq 10;keys:={Ky([Evaluate(PL![Rat(v):v in pt[2]],r[1]),Evaluate(PL![Rat(v):v in pt[3]],r[1])],d):r in roots};
 if idx eq 30 then own:=[r:r in source[2]|r[1] eq 5];sm:=[r:r in source[3]|r[1] eq 5];expected:={own[1][3],sm[19][3]};else own:=[r:r in source[2]|r[1] eq 10];expected:={own[54][3]};end if;assert keys eq expected;
 rs:=[K!0,(13*Q)[1]];Append(~rs,-c-rs[2]);done:={};
 for cert in row[4] do
  i,j,ell,r0,v0:=Explode(cert);assert i ne j and i in [0..2] and j in [0..2];k:=GF(ell);PP:=PolynomialRing(k);rr:=k!r0;fp:=PP!f;assert Evaluate(fp,rr) eq 0 and Degree(GCD(fp,Derivative(fp))) eq 0;v:=Evaluate(PP!Eltseq(rs[j+1]-rs[i+1]),rr);assert v eq k!v0 and v ne 0 and v^((ell-1) div 2) eq -1;Include(~done,[i,j]);
 end for;assert #done eq 6;Append(~out,[*idx,Degree(K),52,d,#keys,#done*]);print "Q226_COMPLETE_NATIVE_CASE",idx,"full_torsion=2x26","class_bound=true","six_nonsquares=true",Cputime(tt);
end for;
Write("work/deg10/C226_resolution_20261030/complete_native_result.m",Sprint(out):Overwrite:=true);print "Q226_COMPLETE_RESOLUTION_CERTIFICATE_DONE",Cputime(tt);quit;
