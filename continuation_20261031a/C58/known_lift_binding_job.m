SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(4*10^9);tt:=Cputime();AttachSpec("work/degree9/shared/mdmagma/v2/mdmagma.spec");
k:=GF(3);
function NewKey(xy)
 L:=Parent(xy[1]);dd:=Degree(L,k);vv:=[xy[1],xy[2],xy[1]+xy[2],xy[1]-xy[2]];
 zz:=[v:v in vv|Degree(MinimalPolynomial(v,k)) eq dd];assert #zz gt 0;th:=zz[1];
 BM:=Matrix(k,[Eltseq(th^i,k):i in [0..dd-1]]);IB:=BM^-1;
 return Sprint([Eltseq(MinimalPolynomial(th,k)),Eltseq(Vector(k,Eltseq(xy[1],k))*IB),Eltseq(Vector(k,Eltseq(xy[2],k))*IB)]);
end function;
load "work/deg10/continuation_20261031a/C58/degree10_census.m";
R<a>:=PolynomialRing(Rationals());models:=[* *];
Append(~models,[*R![-1,0,5,-8,4,-4,7,-5,2,-2,1],R![2697695/97,-7865518/97,9444429/97,-5954840/97,6571299/97,-8368652/97,5515994/97,-2594047/97,2167854/97,-925261/97],R![19044/97,-55602/97,66805/97,-42151/97,46488/97,-59213/97,39068/97,-18361/97,15341/97,-6557/97]*]);
Append(~models,[*R![1,-7,20,-37,51,-56,47,-31,15,-5,1],R![44,-435,1759,-3899,5512,-5368,3787,-1868,616,-112],R![-3,34,-135,291,-405,390,-273,134,-44,8]*]);

rows:=[* *];allindices:={Integers()|};
for z in [1..2] do
 ff,bb,cc:=Explode(models[z]);assert Degree(ff) eq 10 and IsIrreducible(ff);K<aa>:=NumberField(ff);b:=Evaluate(bb,aa);c:=Evaluate(cc,aa);E:=EllipticCurve([1-c,-b,-b,K!0,K!0]);P:=E![0,0];assert 29*P eq E!0 and P ne E!0;
 fp:=Polynomial(GF(3),ff);assert IsIrreducible(fp);L<ar>:=ext<GF(3)|fp>;bv:=Evaluate(bb,ar);cv:=Evaluate(cc,ar);EE:=EllipticCurve([L|1-cv,-bv,-bv,0,0]);PP:=EE![0,0];assert 29*PP eq EE!0 and PP ne EE!0;
 ids:={Integers()|};for a in [1..14] do ky:=NewKey(TateNormalForm_xy(EE,a*PP));i:=Index(PointKeys,ky);assert i gt 0;Include(~ids,i);end for;
 orbitindex:=[i:i in [1..#Orbits]|Seqset(Orbits[i]) eq ids];assert #orbitindex eq 1;rep:=Representatives[orbitindex[1]];assert rep in [462,1748];Include(~allindices,rep);
 odd:=[* *];for q in PrimesUpTo(199) do
  if q eq 2 or &or[Denominator(v) mod q eq 0:v in Coefficients(ff) cat Coefficients(bb) cat Coefficients(cc)] then continue;end if;
  fQ:=Polynomial(GF(q),ff);for rr in Roots(fQ) do
   if rr[2] ne 1 then continue;end if;qq:=rr[1];bq:=Evaluate(bb,qq);cq:=Evaluate(cc,qq);Eq:=EllipticCurve([GF(q)|1-cq,-bq,-bq,0,0]);if Discriminant(Eq) eq 0 then continue;end if;
   card:=#Eq;if IsOdd(card) then odd:=[*q,Integers()!qq,card*];break;end if;
  end for;if #odd gt 0 then break;end if;
 end for;assert #odd gt 0;
 Append(~rows,[*z,rep,Sort(Setseq(ids)),odd*]);
 printf "D31A_C58_KNOWN_LIFT model=%o representative=%o orbit_size=%o odd_reduction=%o time=%o\n",z,rep,#ids,odd,Cputime(tt);
end for;assert allindices eq {462,1748};
Write("work/deg10/continuation_20261031a/C58/known_lift_binding.m",Sprintf("KnownLiftBinding:=%o;",rows));
print "D31A_C58_KNOWN_LIFT_BINDING_COMPLETE two_orbits=true exact29torsion=true no2torsion=true";quit;
