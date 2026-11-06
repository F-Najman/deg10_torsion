SetColumns(0); SetQuitOnError(true); SetNthreads(1); SetMemoryLimit(12*10^9); tt:=Cputime();
F:=GF(2); K<t>:=FunctionField(F); P<v>:=PolynomialRing(K);
jt:=(t^2+5*t+13)*(t^4+7*t^3+20*t^2+19*t+1)^3/t;
f:=v^8+jt*v^7+v^4+1; assert IsIrreducible(f); AF<vv>:=FunctionField(f);
assert Genus(AF) eq 7; jaf:=AF!jt;
cs:=Poles(jaf); widths:=[-Valuation(jaf,c):c in cs];
assert #cs eq 4 and &and[Degree(c) eq 1:c in cs] and Seqset(widths) eq {1,7,13,91};
J:=JZero(91); fj:=FrobeniusPolynomial(J,2); if Type(fj) eq SeqEnum then fj:=&*fj; end if; nj:=Integers()!Evaluate(fj,1);
printf "C91_BASE_READY genus=%o widths=%o Jorder=%o time=%o\n",Genus(AF),widths,nj,Cputime(tt);
DG:=DivisorGroup(AF);Dzero:=DG!0;G,gm:=RayClassGroup(Dzero);
inv:=Invariants(G);assert &*[a:a in inv|a ne 0] eq nj and #[a:a in inv|a eq 0] eq 1;
printf "C91_PICARD invariants=%o time=%o\n",inv,Cputime(tt);
cols:=[];
for i in [1..4] do
 Append(~cols,Eltseq(Divisor(cs[i])@@gm));
 printf "C91_CUSP_CLASS index=%o class=%o time=%o\n",i,cols[i],Cputime(tt);
end for;
expect:=[w in {91,13} select 0 else 3:w in widths];
chars:=[];
for r in CartesianPower({0..3},Ngens(G)) do
 rr:=[r[i]:i in [1..Ngens(G)]];
 if not &and[inv[i]*rr[i] mod 4 eq 0:i in [1..#inv]] then continue;end if;
 if &and[(&+[cols[i][j]*rr[j]:j in [1..#inv]]) mod 4 eq expect[i]:i in [1..4]] then Append(~chars,rr);end if;
end for;
printf "C91_C4_CUSP_CHARACTERS count=%o images=%o time=%o\n",#chars,chars,Cputime(tt);

print "C91_CHARACTER_PROBE_COMPLETE"; quit;
