SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(6*10^9);tt:=Cputime();
root:="work/deg10/next_cases_20261030/C242/direct_hecke17/";input:=eval Read(root cat INPUT_FILE);fields:=input[1];cases:=input[2];assert ELL in [11,17] and SHARD ge 0 and SHARD lt SHARDS;
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");load "work/deg10/next_cases_20261030/C242/direct_hecke17/deps/pair_af_place_linear.m";
Y:=MDX1(42,GF(5));C:=Curve(Y);assert Genus(C) eq 25;cs:=Cusps(Y);os:=CuspOrbitsQ(Y);c0:=Support([o:o in os|Degree(o) eq 1][1])[1];AF,am:=AlgorithmicFunctionField(FunctionField(C));coords:=[am(FunctionField(C).i):i in [1,2]];
function IsogeniesProof(E, points, q)
    assert IsPrime(q) and IsOdd(q) and q ne Characteristic(BaseRing(E));
    K := BaseRing(E);
    fq := DivisionPolynomial(E, q);
    fac := [t[1] : t in Factorization(fq)];
    assert &+[Degree(h) : h in fac] eq (q^2-1) div 2;   // separable: distinct x-coordinates
    a1, a2, a3, a4, a6 := Explode(aInvariants(E));
    unassigned := [true : h in fac];
    result := [* *];
    total := 0;
    for i in [1..#fac] do
        if not unassigned[i] then continue; end if;
        h := fac[i]; d := Degree(h);
        if d eq 1 then F1 := K; x1 := -Coefficient(h, 0)/Coefficient(h, 1);
        else F1 := ext<K | h>; x1 := F1.1; end if;
        // y-coordinate: solve y^2 + a1 x y + a3 y = x^3 + a2 x^2 + a4 x + a6 over F1 or a quadratic extension
        PY<Y> := PolynomialRing(F1);
        ypol := Y^2 + (F1!a1*x1 + F1!a3)*Y - (x1^3 + F1!a2*x1^2 + F1!a4*x1 + F1!a6);
        rts := Roots(ypol);
        if #rts gt 0 then F2 := F1; y1 := rts[1][1];
        else F2 := ext<F1 | 2>; y1 := Roots(PolynomialRing(F2)!ypol)[1][1]; end if;
        E2 := EllipticCurve([F2!a1, F2!a2, F2!a3, F2!a4, F2!a6]);
        P1 := E2![F2!x1, y1];
        assert q*P1 eq E2!0 and P1 ne E2!0;
        xs := [ F1!((k*P1)[1]) : k in [1..(q-1) div 2] ];          // x-coordinates lie in F1
        PX<X> := PolynomialRing(F1);
        g := &*[X - xk : xk in xs];                                  // kernel polynomial of G = <P1>
        cs := Coefficients(g);
        L := sub< F1 | cs cat [F1!K.1] >;                            // field of definition of G, contains K
        e := Degree(L) div Degree(K);
        // mark the factors belonging to this Frobenius-orbit: divisors of Norm_{L/K}(g)
        qK := #K;
        Ng := PX!1; gg := g;
        for j in [1..e] do
            Ng *:= gg;
            gg := PX![c^qK : c in Coefficients(gg)];
        end for;
        assert gg eq g;   // Frobenius^e fixes g
        NgK := PolynomialRing(K)![K!c : c in Coefficients(Ng)];
        deg_assigned := 0;
        for j in [1..#fac] do
            if unassigned[j] and IsDivisibleBy(NgK, fac[j]) then unassigned[j] := false; deg_assigned +:= Degree(fac[j]); end if;
        end for;
        assert deg_assigned eq e*(q-1) div 2;
        // the isogeny over L
        PL<XL> := PolynomialRing(L);
        gL := PL![L!c : c in cs];
        EL := EllipticCurve([L!(F1!a) : a in aInvariants(E)]);
        EG, phi := IsogenyFromKernel(EL, gL);
        images := [phi(EL![L!(F1!(K!P[1])), L!(F1!(K!P[2]))]) : P in points];
        Append(~result, <EG, images, e, gL, EL, [EL![L!(F1!(K!P[1])), L!(F1!(K!P[2]))]:P in points], L!(F1!K.1)>);
        total +:= e;
    end for;
    assert total eq q+1;
    return result;
end function;
function FieldData(K) return [Integers()!a:a in Coefficients(DefiningPolynomial(K,GF(5)))];end function;
function ElementData(v) return [Integers()!a:a in Eltseq(v,GF(5))];end function;
function CurveData(E,P,Q) return [*FieldData(BaseRing(E)),[ElementData(a):a in aInvariants(E)],[ElementData(a):a in Eltseq(P)],[ElementData(a):a in Eltseq(Q)]*];end function;
function CoeffData(f) return [[Integers()!a:a in Coefficients(Numerator(f))],[Integers()!a:a in Coefficients(Denominator(f))]];end function;
function FData(f) return [CoeffData(c):c in Eltseq(AF!f)];end function;
OF:=MaximalOrderFinite(AF);OI:=MaximalOrderInfinite(AF);bf:=[AF!b:b in Basis(OF)];bi:=[AF!b:b in Basis(OI)];
Write(root cat TAG cat "/model_bases.m",Sprint([*Degree(AF),[CoeffData(c):c in Coefficients(DefiningPolynomial(AF))],[FData(c):c in coords],[FData(b):b in bf],[FData(b):b in bi]*]):Overwrite:=true);
R0:=PolynomialRing(GF(5));fs:=AssociativeArray();for fd in fields do fs[fd[1]]:=ext<GF(5)|R0!fd[2]>;end for;
cache:=AssociativeArray();pointProof:=[* *];records:=[* *];remaining:=[* *];
procedure APoint(~answer,~cache,~pointProof,d,ky,mapno)
 kk:=[mapno,d] cat ky;if IsDefined(cache,kk) then answer:=cache[kk];return;end if;
 K:=fs[d];z:=K.1;b:=Evaluate(R0!ky[1..d],z);c:=Evaluate(R0!ky[d+1..2*d],z);E:=EllipticCurve([K|0,c,0,b*(1-b-c),0]);P:=E![0,0];Q:=E![b,b];assert 42*Q eq E!0 and 21*Q ne E!0 and 14*Q ne E!0 and 6*Q ne E!0 and P ne 21*Q;
 if mapno eq 1 then EG:=E;QG:=Q;elif mapno eq 2 then EG:=E;QG:=P+Q;elif mapno eq 3 then EG:=E;QG:=P+22*Q;else T:=mapno eq 4 select P else mapno eq 5 select 21*Q else P+21*Q;U:=mapno eq 5 select P+22*Q else Q;EG,phi2:=IsogenyFromKernel(E,PolynomialRing(K).1-T[1]);QG:=phi2(U);end if;
 assert 42*QG eq EG!0 and 21*QG ne EG!0 and 14*QG ne EG!0 and 6*QG ne EG!0;
 y:=C48PairAFPlaceLinear(coords,TateNormalForm_xy(EG,QG));yd:=C48PairAFPlaceLinear(coords,TateNormalForm_xy(EG,ELL*QG));de:=Degree(y);assert de eq Degree(yd) and d mod de eq 0;weight:=d div de;
 isog:=IsogeniesProof(EG,[QG,QG],ELL);assert &+[q[3]:q in isog] eq ELL+1;Dz:=DivisorGroup(AF)!0;idata:=[* *];
 for q in isog do
  EH,pts,e,gg,Esource,pointsource,embedK:=Explode(q);vv:=TateNormalForm_xy(EH,pts[1]);yy:=C48PairAFPlaceLinear(coords,vv);dd:=Degree(yy);assert d*e mod dd eq 0;mult:=d*e div dd;Dz+:=mult*Divisor(yy);
  Append(~idata,[*e,CurveData(Esource,pointsource[1],pointsource[2]),[ElementData(a):a in Coefficients(gg)],CurveData(EH,pts[1],pts[2]),ElementData(embedK),[ElementData(a):a in vv],dd,mult*]);
 end for;
 assert Degree(Dz) eq (ELL+1)*d;answer:=Dz-weight*Divisor(y)-ELL*weight*Divisor(yd);assert Degree(answer) eq 0;cache[kk]:=answer;
 Append(~pointProof,[*kk,CurveData(E,P,Q),CurveData(EG,QG,QG),de,weight,idata*]);Write(root cat TAG cat "/explicit_hecke_data.m",Sprint(pointProof):Overwrite:=true);print "C242_DIRECT_HECKE_POINT_READY",ELL,TAG,mapno,d,#pointProof,Cputime(tt);
end procedure;
for localid in [1..#cases] do
 if (localid-1) mod SHARDS ne SHARD then continue;end if;ci:=cases[localid][1];terms:=cases[localid][2];cd:=10-&+[t[1]:t in terms];assert cd ge 0 and cd le 6;
 if ELL eq 11 and cd ne 0 then Append(~records,[*localid,ci,cd,"SKIPPED_NONZERO_CUSP_COMPLEMENT",[* *],true*]);Append(~remaining,cases[localid]);continue;end if;
 passed:=true;checks:=[* *];
 for mapno in [1..6] do
  Adiv:=DivisorGroup(AF)!0;for term in terms do aa:=DivisorGroup(AF)!0;APoint(~aa,~cache,~pointProof,term[1],term[2],mapno);Adiv+:=aa;end for;assert Degree(Adiv) eq 0;
  print "C242_DIRECT_HECKE_RR_START",ELL,TAG,ci,mapno,Cputime(tt);bbasis:=Basis(Adiv:Reduction:=false,Simplification:="None");assert #bbasis in [0,1];if #bbasis eq 1 then assert Divisor(bbasis[1])+Adiv eq DivisorGroup(AF)!0;end if;
  I:=ideal<OF|1>;J:=ideal<OI|1>;ss,mm:=Support(Adiv);for k in [1..#ss] do if IsFinite(ss[k]) then I*:=Ideal(ss[k])^(-mm[k]);else J*:=Ideal(ss[k])^(-mm[k]);end if;end for;
  fname:=Sprintf("ideal_basis_case%o_map%o.m",ci,mapno);Write(root cat TAG cat "/" cat fname,Sprint([*ci,mapno,[FData(b):b in Basis(I)],[FData(b):b in Basis(J)]*]):Overwrite:=true);
  Append(~checks,[*mapno,#bbasis,fname,Sprint(Adiv),[FData(f):f in bbasis]*]);print "C242_DIRECT_HECKE_CASE_MAP_DONE",ELL,TAG,ci,mapno,#bbasis,Cputime(tt);
  if #bbasis eq 0 then passed:=false;break;end if;
 end for;
 Append(~records,[*localid,ci,cd,"TESTED",checks,passed*]);if passed then Append(~remaining,cases[localid]);end if;
 Write(root cat TAG cat "/arithmetic.m",Sprint([*ELL,fields,records,remaining*]):Overwrite:=true);
end for;
Write(root cat TAG cat "/arithmetic.m",Sprint([*ELL,fields,records,remaining*]):Overwrite:=true);print "C242_DIRECT_HECKE_COMPLETE",ELL,TAG,#records,#remaining,Cputime(tt);quit;
