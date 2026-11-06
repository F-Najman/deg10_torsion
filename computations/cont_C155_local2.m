SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(2*10^9);
assert 16 mod 5 eq 1 and 16 mod 31 eq 16 and 16^5 mod 155 eq 1;
assert &and[16^i mod 155 notin [1,154]:i in [1..4]];
H:=Sort(Setseq({r mod 31:r in [1,2,4,8,16,-1,-2,-4,-8,-16]}));assert #H eq 10;
idx:=EulerPhi(155)*155*&*[Rationals()|1+1/p[1]:p in Factorization(155)]/2;
assert idx eq 11520 and 80 lt 325*idx/32768;
A:=JH(31,[2]);assert Dimension(A) eq 6;
load "work/deg10/inputs/dcm_list_all_min_orders_X1_10k.m";assert dcm_list_all_min_orders_X1[155][2] eq 120;
rows:=[];
for q in [2,4] do
 F:=GF(q);counts:={Integers()|};models:=0;
 for a1,a2,a3,a4,a6 in F do
  try E:=EllipticCurve([a1,a2,a3,a4,a6]);catch err continue;end try;
  s:=#E;if s mod 5 ne 0 then continue;end if;
  models+:=1;t:=q+1-s;Include(~counts,t);
  assert s eq 5 and not IsSquare(GF(31)!(t^2-4*q));
 end for;
 assert #counts eq 1;Append(~rows,[q,Rep(counts),models]);print "DECIC155_LOCAL_TWO_NO_LINE",q,Rep(counts),models;
end for;
Write("work/deg10/continuation_20261029/C155_local2_arithmetic.txt",Sprint([*H,idx,325*idx/32768,120,rows*]));
print "DECIC155_LOCAL_TWO_COMPLETE";quit;
